# Define structure
mutable struct inputcom
    y        :: Vector{Float64}
    x        :: Vector{Float64}
    w        :: Array{Float64, 2}
    H_min    :: Int64
    H_max    :: Int64
    type     :: String
    r        :: Float64
    λ        :: Float64
  end

mutable struct outputcom
    T      :: Int64
    H_min  :: Int64
    H_max  :: Int64
    HR     :: Int64
    K      :: Int64
    B      :: Array{Float64, 2}
    P      :: Array{Float64, 2}
    λ      :: Float64
    type   :: String
    δ      :: Float64
    idx    :: Array{Float64, 2}
    Θ      :: Vector{Float64}
    IR     :: Array{Float64, 2}
    X      :: Array{Float64, 2}
    Y      :: Vector{Float64}
    isreg  :: Bool
end

#Define functional form
function initalz(df, indx, indx_str, param)
    T, k = size(df)
    ind_response, ind_shock = indx[1], indx[2];
    ind_all = 1:k;

    H_min, H_max = 1, indx[3];
    P = indx[4]
    r, λ = param[1], param[2];

    # Indicators for variables
    ind_contempo = filter(x -> x ≠ ind_shock, ind_all);

    y  = df[:, ind_response]; # endogenous variable
    x  = df[:, ind_shock]; # endoegnous variable related to the shock

    # control variables (contemporaneous vars, lagged vars)
    w  = [ df[:, ind_contempo]  matlag( df , P ) ];
    w[.!isfinite.(w)] .= 0;

    return inputcom(y, x, w, H_min, H_max, indx_str, r, λ);
end

function slp(packagedinput)
    # Unpackage data
    y        = packagedinput.y;
    x        = packagedinput.x;
    w        = packagedinput.w;
    H_min    = packagedinput.H_min;
    H_max    = packagedinput.H_max;
    type     = packagedinput.type;
    r        = packagedinput.r;
    λ        = packagedinput.λ;

    # Assign δ
    if isempty(w)
        δ = std(x);
    else
        δ = std(x - w * inv(w' * w) * w' * x);
    end

    # Is it the regular regression-type?
    isreg = type == "reg";

    T  = length(y);
    HR = H_max + 1 - H_min;

    # construct the B-spline basis functions
    κ = 3;
    if ~isreg
        B = bspline((H_min:H_max)', H_min, (H_max + 1), (H_max + 1 - H_min), κ);
        K = size(B, 2);
    else
        K = HR;
        B = fill(0.0, (H_max, H_max + κ));
        λ = 0;
    end

    # building up the regression representation of the local projection
    idx = nan((H_max+1)*T, 2);
    Y   = nan((H_max+1)*T, 1);
    Xb  = zeros(Float64, (H_max+1)*T, K);
    Xc  = zeros(Float64, (H_max+1)*T, HR, size(w,2)+1 );
    # NOTE: For potential for parallel computing, which may require to donwsize
    #   bit size, will lead to substantially inaccurate results for inverse functions
    # https://stackoverflow.com/questions/51562302/inverse-matrix-results-different-in-matlab-and-python
    w = hcat(ones(T,1), w);

    for t ∈ 1:(T - H_min)

        idx_beg = (t - 1)*HR + 1;
        idx_end = t * HR;

        idx[idx_beg:idx_end, 1] .= t;
        idx[idx_beg:idx_end, 2] .= H_min:H_max;

        # y
        y_range = (t + H_min) : min((t + H_max), T)';
        Y[idx_beg:idx_end] = [ y[y_range] ; nan(HR-length(y_range),1) ];

        # x
        if isreg
            Xb[ idx_beg:idx_end , : ] = eye(HR) .* x[t];
        else
            Xb[ idx_beg:idx_end , : ] = B * x[t];
        end

        # w
        for i ∈ 1:size(w,2)
            Xc[idx_beg:idx_end , : , i ] = eye(HR)*w[t,i];
        end
    end

    X = Xb;
    for i ∈ 1:size(w,2)
        X = hcat(X, Xc[:, :, i]);
    end

    select = isfinite.(Y)[:];
    idx = idx[select, :];
    Y   = Y[select];
    X   = X[select,:];

    IR = zeros(Float16, H_max+1, 1);

    if isreg
        Θ = (X' * X)\X'*Y
        IR[(H_min+1):end] = Θ[1:K] .* δ;

        # Parameter not used but still defined,
        P = zeros(Float16, size(X,2), size(X,2) );
    else
        P = zeros(Float16, size(X,2), size(X,2));

        D = eye(K);
        for k = 1:r #k ∈ 1:r
            D = diff(D, dims = 1);
        end

        P[1:K,1:K] = D' * D;

        Θ = ( (X' * X) + (λ .* P) )\( X'*Y );

        IR[(1+H_min):end] = B * Θ[1:K] .* δ;
    end

    output = outputcom(T, H_min, H_max, HR, K, B,P,λ,type,δ,idx,Θ,IR,X,Y,isreg);

    return(output);
    # debug stuff
    # Bs is for display / debugging pourposes only
    # Bs = bspline( (H_min:0.1:H)' , H_min , H+1 , H+1-H_min , 3 );
    # obj.Bs = Bs;
end

function slp_ci(outputcomⱼ)
    # Confidence band imputation, input from outcome
    X = outputcomⱼ.X
    Y = outputcomⱼ.Y
    
    
    # Parameter saved
    λ = outputcomⱼ.λ
    P = outputcomⱼ.P

    H_max = outputcomⱼ.H_max;
    H_min = outputcomⱼ.H_min;
    T = outputcomⱼ.T;
    n_param = length(outputcomⱼ.Θ);
    ω = vcat(0.5, ((H_max + 1) .- (1:H_max))./(H_max + 1));
    idx = outputcomⱼ.idx;
    K = outputcomⱼ.K
    δ = outputcomⱼ.δ
    
    B = outputcomⱼ.B
    isreg = outputcomⱼ.isreg;

    # Point estimate
    XXₚ = ((X' * X) + (λ .* P));
    Θ = outputcomⱼ.Θ;
    U = Y .- X * Θ;
    
    HR = outputcomⱼ.HR;
    
    V = zeros(n_param, n_param);
    titer = 0;
    S1 = X[ idx[:,1] .== titer , : ]' * U[ idx[:,1] .== titer ]
    S2 = X[ idx[:,1] .== (titer - 1) , : ]' * U[ idx[:,1] .== (titer - 1) ]
    for iter ∈ 0:H_max
        ggt = zeros(n_param, n_param);
        for titer ∈ (iter + 1):(T - HR - 1)
            S1 = X[ idx[:,1] .== titer , : ]' * U[ idx[:,1] .== titer ];
            S2 = X[ idx[:,1] .== (titer - iter) , : ]' * U[ idx[:,1] .== (titer - iter) ];
            ggt += S1 * S2' + S2 * S1';
        end
        V += ω[iter + 1] .* ggt;
    end

    VC = XXₚ^(-1) * V * XXₚ^(-1);
    
    # To save confidence band
    conf = nan(H_max + 1, 2);

    # 90% confidence band;
    t_lb = quantile(Normal(0,1), 0.05)
    t_ub = quantile(Normal(0,1), 0.95)
    trail_Θ = Θ[1:K];

    if isreg
        h = 1:(H_max + 1 - H_min);
        ρ = sqrt.(diag(VC[h,h]));
        conf[(1+H_min):end, 1] = (trail_Θ .* δ) .+ (ρ .* (δ * t_lb))
        conf[(1+H_min):end, 2] = (trail_Θ .* δ ).+ (ρ .* (δ * t_ub))
    else # For spline regression
        h = 1:K;
        ρ = sqrt.(diag(B * VC[h,h] * B'));
        conf[(1+H_min):end, 1] = (B * trail_Θ .* δ) .+ (ρ .* (δ * t_lb))
        conf[(1+H_min):end, 2] = (B * trail_Θ .* δ) .+ (ρ .* (δ * t_ub))
    end
    # Initialize the first horizon - impact
    conf[1, 1] = outputcomⱼ.IR[1];
    conf[1, 2] = outputcomⱼ.IR[1];
    return(conf);
end

function slpᵥ(output, λₘ = 10)
    X = output.X;
    Y = output.Y;
    P = output.P;
    T = output.T;

    λᵥ = Array(1.0 : 0.05 : λₘ) .* T;

    L = length(λᵥ);
    resvec = zeros(Float64, L, 2);
    resvec[:, 1] = λᵥ[:];

    for l ∈ 1:L
        S = X * inv((X' * X) + (λᵥ[l] .* P)) * X';
        #X * (((X' * X) + (λᵥ[l] .* P)) \ X')
        # Need to convert
        resvec[l,2] = rss(Y, S*Y, (1 .- diag(S)));
    end

    λₒ = resvec[(minimum(resvec[:,2]) .== resvec[:,2]), 1][];

    return(λₒ, resvec);
end

function bspline(x, xl, xr, ndx, bdeg)

    dx = (xr - xl) / ndx;

    t = xl .+ dx .* vec(-bdeg:(ndx - 1))';
    T = ones(Float16, length(x), 1) * t;
    X = x' * ones(Float16, 1, length(t));
    P = (X .- T) ./ dx;

    B = (T .<= X) .& (X .< (T .+ dx));
    r = vcat(2:length(t), 1);

    for k ∈ 1:bdeg
        B = (P .* B + (k + 1 .- P) .* B[:, r]) ./ k;
    end
    return B;
end

function rss(y, ŷ, σ)
    res = sum(((y .- ŷ) ./ σ).^2);
    return res;
end

function matlag(x, lag)
    t, k = size(x);
    res = zeros(Float16, (t-lag), lag * k);

    for i ∈ 1:lag
        res[:, ((i-1)*k + 1):(i*k)] = x[(lag-i + 1):end-i ,:];
    end

    res = nan(t, lag * k);
    for i ∈ 1:lag
        lag_array = x[1:(end-i),:];
        res[(i+1):end, ((i-1)*k + 1):(i*k)] = lag_array;
    end

    return res;
end

eye(n) = Matrix{Float16}(I, n, n);

function nan(x, y)
    temp = zeros(Float16, x, y);
    return replace(temp, 0 => NaN);
end