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
        δ = std( x .- w*inv(w'*w)*w'*x );
    end

    # Is it the regular regression-type?
    isreg = type == "reg";


    T  = length(y);
    HR = H_max + 1 - H_min;

    # construct the B-spline basis functions
    if ~isreg
        B = bspline( (H_min:H_max)' , H_min , H_max+1 , H_max+1-H_min , 3 );
        K = size( B , 2 );
    else
        K = HR;
    end

    # building up the regression representation of the local projection
    idx = nan( (H_max+1)*T , 2 );
    Y   = nan( (H_max+1)*T , 1 );
    Xb  = zeros( (H_max+1)*T , K );
    Xc  = zeros( (H_max+1)*T , HR , size(w,2)+1 );

    w = hcat(ones(T,1), w);

    for t = 1:(T - H_min)

        idx_beg = (t - 1)*HR + 1;
        idx_end = t * HR;

        idx[idx_beg:idx_end, 1] = t;
        idx[idx_beg:idx_end, 2] = H_min:H_max;

        # y
        y_range = (t + H_min) : min((t + H_max), T)';
        Y[idx_beg:idx_end] = [ y[y_range] ; nan(HR-length(y_range),1) ];

        # x
        if isreg
            Xb[ idx_beg:idx_end , : ] = eye(HR)*x[t];
        else
            Xb[ idx_beg:idx_end , : ] = B*x[t];
        end

        # w
        for i = 1:size(w,2)
            Xc[idx_beg:idx_end , : , i ] = eye(HR)*w[t,i];
        end

    end


    # STOP line 59 of matlab
end

function bspline(x, xl, xr, ndx, bdeg)
    dx = (xr - xl) / ndx;

    t = xl .+ dx .* [-bdeg:(ndex - 1)];
    T = (0 .*x + 1) .* t;

    X = x * (0 .* t .+ 1);
    P = (X - T) / dx;

    B = (T <= X) & (X < (T + dx));
    r = [2:length(t) 1];

    for k = 1:bdeg
        B = (P .* B + (k + 1 - P) .* B(:, r)) / k;
    end
    return B;
end

df_test = randn(20,7)
function matlag(x, lag)
    t, k = size(x);
    res = zeros((t-lag), lag * k)
    for i = 1:lag
        res[:, ((i-1)*k + 1):(i*k)] = x[(lag-i + 1):end-i ,:]
    end
    return res
end

eye(n) = Matrix{Float64}(I, n, n)

function nan(x, y)
    temp = zeros(x, y);
    return replace(temp, 0 => NaN)
end
