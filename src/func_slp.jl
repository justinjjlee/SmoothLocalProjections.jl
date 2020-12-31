mutable struct outputcom
    T
    H_min
    H_max
    HR
    K
    B
    P
    λ
    type
    δ
    idx
    Θ
    IR
end

function smoothlocalprojection(packagedinput)
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
    Xb  = zeros(Float16, (H_max+1)*T, K);
    Xc  = zeros(Float16, (H_max+1)*T, HR, size(w,2)+1 );

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
        B = 0;
        P = zeros(Float16, size(X,2), size(X,2) );
        λ = 0;
    else
        P = zeros(Float16, size(X,2), size(X,2));

        D = eye(K);
        for k ∈ 1:r
            D = diff(D, dims = 1);
        end

        P[1:K,1:K] = D' * D;

        Θ = ( (X' * X) + (λ .* P) )\( X'*Y );

        IR[(1+H_min):end] = B * Θ[1:K] .* δ;
    end

    output = outputcom(T, H_min, H_max, HR, K, B,P,λ,type,δ,idx,Θ,IR) ;

    return(output);
    # debug stuff
    # Bs is for display / debugging pourposes only
    # Bs = bspline( (H_min:0.1:H)' , H_min , H+1 , H+1-H_min , 3 );
    # obj.Bs = Bs;
end

function bspline(x, xl, xr, ndx, bdeg)

    dx = (xr - xl) / ndx;

    t = xl .+ dx .* vec(-bdeg:(ndx - 1))'
    T = ones(Float16, length(x), 1) * t;
    X = x' * ones(Float16, 1, length(t));
    P = (X .- T) ./ dx;

    B = (T .<= X) .& (X .< (T .+ dx));
    r = vcat(2:length(t), 1);

    for k = 1:bdeg
        B = (P .* B + (k + 1 .- P) .* B[:, r]) ./ k;
    end
    return B;
end

function matlag(x, lag)
    t, k = size(x);
    res = zeros(Float16, (t-lag), lag * k)

    for i = 1:lag
        res[:, ((i-1)*k + 1):(i*k)] = x[(lag-i + 1):end-i ,:]
    end

    res = nan(t, lag * k);
    for i = 1:lag
        lag_array = x[1:(end-i),:];
        res[(i+1):end, ((i-1)*k + 1):(i*k)] = lag_array;
    end

    return res
end

eye(n) = Matrix{Float16}(I, n, n)

function nan(x, y)
    temp = zeros(Float16, x, y);
    return replace(temp, 0 => NaN)
end
