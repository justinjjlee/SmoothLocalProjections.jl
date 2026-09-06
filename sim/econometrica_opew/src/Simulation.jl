module Simulation

using LinearAlgebra
using Statistics
using Distributions

export VARModel, estimate_var, select_lag_aic, estimate_lp, simulate_data, compute_population_irfs

struct VARModel
    p::Int
    n::Int
    B::Matrix{Float64}         # n × (n*p + 1) [slopes, intercept]
    A_comp::Matrix{Float64}    # np × np companion matrix
    Sigma_u::Matrix{Float64}   # n × n residual covariance
    C_chol::Matrix{Float64}    # n × n lower Cholesky factor
    nu::Vector{Float64}        # n × 1 normalized impact vector for shock 1
    nu_comp::Vector{Float64}   # np × 1 normalized impact vector
    resid::Matrix{Float64}     # (T - p) × n residuals
    X_expand::Matrix{Float64}  # (T - p) × (np + 1) regressor matrix
end

"""
    estimate_var(Y::Matrix{Float64}, p::Int)

Estimate reduced-form VAR(p) with intercept via OLS.
"""
function estimate_var(Y::Matrix{Float64}, p::Int)
    T_raw, n = size(Y)
    T_eff = T_raw - p
    
    Y_lhs = Y[(p + 1):end, :]
    X = zeros(T_eff, n * p)
    for l in 1:p
        X[:, ((l - 1)*n + 1):(l*n)] = Y[(p - l + 1):(end - l), :]
    end
    X_expand = hcat(X, ones(T_eff, 1))
    
    betahat = Matrix((X_expand \ Y_lhs)') # n × (n*p + 1)
    resid = Y_lhs - X_expand * betahat'
    Sigma_u = (resid' * resid) ./ (T_eff - size(betahat, 2))
    
    # Companion matrix
    np = n * p
    A_comp = zeros(np, np)
    A_comp[1:n, :] = betahat[:, 1:end-1]
    if p > 1
        A_comp[(n + 1):end, 1:(n*(p - 1))] = Matrix(1.0I, n*(p - 1), n*(p - 1))
    end
    
    # Cholesky factor with lower triangular normalization
    C_chol = cholesky(Hermitian(Sigma_u)).L
    nu = C_chol[:, 1] ./ C_chol[1, 1]
    nu_comp = zeros(np)
    nu_comp[1:n] = nu
    
    return VARModel(p, n, betahat, A_comp, Sigma_u, C_chol, nu, nu_comp, resid, X_expand)
end

"""
    compute_var_irf_and_se(model::VARModel, H_max::Int=50, i_star::Int=7)

Compute VAR impulse response of variable i_star to shock 1, and
homoskedastic delta-method standard errors.
"""
function compute_var_irf_and_se(model::VARModel, H_max::Int=50, i_star::Int=7)
    np = size(model.A_comp, 1)
    n = model.n
    p = model.p
    T_eff = size(model.X_expand, 1)
    
    e_i = zeros(np); e_i[i_star] = 1.0
    irfs = zeros(H_max + 1)
    ses = zeros(H_max + 1)
    
    # Precompute powers of A_comp
    A_pows = [Matrix(1.0I, np, np)]
    for h in 1:H_max
        push!(A_pows, A_pows[end] * model.A_comp)
    end
    
    # V_A = inv(X'X) ⊗ Sigma_u for companion slope coefficients
    X_slopes = model.X_expand[:, 1:end-1]
    inv_XX = inv(X_slopes' * X_slopes)
    
    for h in 0:H_max
        irfs[h + 1] = dot(e_i, A_pows[h + 1] * model.nu_comp)
        
        if h == 0
            # Horizon 0 response is normalized to 1 for proxy, or impact on variable i
            ses[h + 1] = sqrt(model.Sigma_u[i_star, i_star] / (T_eff * model.C_chol[1, 1]^2))
        else
            # Delta method derivative wrt companion matrix slopes
            # d(e_i' A^h nu) / d(vec(A_1..A_p))
            # Psi_h = ∑_{l=1}^h A^{h-l} nu e_i' A^{l-1}
            Psi_h = zeros(np, np)
            for l in 1:h
                Psi_h .+= A_pows[h - l + 1] * model.nu_comp * (e_i' * A_pows[l])
            end
            # Restrict to first n columns (the derivative with respect to companion slopes B)
            G = Matrix(Psi_h[:, 1:n]') # n × np
            var_h = tr(G * inv_XX * G' * model.Sigma_u)
            ses[h + 1] = sqrt(max(1e-12, var_h))
        end
    end
    
    return irfs, ses
end

"""
    select_lag_aic(Y::Matrix{Float64}, p_max::Int=24)

Select optimal lag length p ∈ 1:p_max minimizing AIC:
AIC(p) = ln|Sigma_p| + 2*p*n^2 / T.
"""
function select_lag_aic(Y::Matrix{Float64}, p_max::Int=24)
    T_raw, n = size(Y)
    best_aic = Inf
    best_p = 1
    
    # Use fixed sample for fair comparison or rolling
    for p in 1:p_max
        T_eff = T_raw - p
        Y_lhs = Y[(p + 1):end, :]
        X = zeros(T_eff, n * p)
        for l in 1:p
            X[:, ((l - 1)*n + 1):(l*n)] = Y[(p - l + 1):(end - l), :]
        end
        X_expand = hcat(X, ones(T_eff, 1))
        betahat = Matrix((X_expand \ Y_lhs)')
        resid = Y_lhs - X_expand * betahat'
        Sigma_p = (resid' * resid) ./ T_eff
        
        # log determinant
        log_det = logdet(Hermitian(Sigma_p))
        aic = log_det + (2.0 * p * (n^2)) / T_raw
        if aic < best_aic
            best_aic = aic
            best_p = p
        end
    end
    return best_p
end

"""
    estimate_lp(Y::Matrix{Float64}, p::Int, H_max::Int=50, i_star::Int=7, j_star::Int=1)

Estimate Local Projection impulse responses and homoskedastic standard errors:
y_{i*, t+h} = β_h y_{j*, t} + γ_h' [y_{t-1}, ..., y_{t-p}] + const + error.
"""
function estimate_lp(Y::Matrix{Float64}, p::Int, H_max::Int=50, i_star::Int=7, j_star::Int=1)
    T_raw, n = size(Y)
    irfs = zeros(H_max + 1)
    ses = zeros(H_max + 1)
    
    for h in 0:H_max
        T_h = T_raw - p - h
        if T_h <= n * p + 2
            break
        end
        
        # LHS: y_{i*, t+h} for t = (p + 1):(T_raw - h)
        y_lhs = Y[(p + 1 + h):(T_raw), i_star]
        
        # Regressor of interest: y_{j*, t}
        shock_reg = Y[(p + 1):(T_raw - h), j_star]
        
        # Lagged controls: y_{t-1}, ..., y_{t-p}
        W_lags = zeros(T_h, n * p)
        for l in 1:p
            W_lags[:, ((l - 1)*n + 1):(l*n)] = Y[(p - l + 1):(T_raw - h - l), :]
        end
        
        # Controls ordered before j* (none if j*=1)
        if j_star > 1
            W_prior = Y[(p + 1):(T_raw - h), 1:(j_star - 1)]
            X_all = hcat(shock_reg, W_prior, W_lags, ones(T_h, 1))
        else
            X_all = hcat(shock_reg, W_lags, ones(T_h, 1))
        end
        
        # OLS
        XtX = X_all' * X_all
        Xty = X_all' * y_lhs
        inv_XtX = inv(XtX)
        beta_all = inv_XtX * Xty
        irfs[h + 1] = beta_all[1]
        
        resid_h = y_lhs - X_all * beta_all
        sigma2_h = sum(resid_h.^2) / (T_h - size(X_all, 2))
        ses[h + 1] = sqrt(max(1e-12, sigma2_h * inv_XtX[1, 1]))
    end
    
    return irfs, ses
end

"""
    compute_population_irfs(pop_model::VARModel, H_max::Int=50, i_star::Int=7)

Compute the exact population impulse response of variable i_star to shock 1.
"""
function compute_population_irfs(pop_model::VARModel, H_max::Int=50, i_star::Int=7)
    np = size(pop_model.A_comp, 1)
    e_i = zeros(np); e_i[i_star] = 1.0
    
    true_irfs = zeros(H_max + 1)
    A_pow = Matrix(1.0I, np, np)
    for h in 0:H_max
        true_irfs[h + 1] = dot(e_i, A_pow * pop_model.nu_comp)
        A_pow = A_pow * pop_model.A_comp
    end
    return true_irfs
end

"""
    simulate_data(pop_model::VARModel, T_sim::Int=720, burn_in::Int=100)

Simulate synthetic time series from VAR(18) with Gaussian shocks.
"""
function simulate_data(pop_model::VARModel, T_sim::Int=720, burn_in::Int=100)
    n = pop_model.n
    p = pop_model.p
    T_tot = T_sim + burn_in
    
    # Innovations
    shocks = randn(T_tot, n) * pop_model.C_chol'
    
    # Simulate
    Y_sim = zeros(T_tot, n)
    for t in (p + 1):T_tot
        y_val = zeros(n)
        for l in 1:p
            y_val .+= pop_model.B[:, ((l - 1)*n + 1):(l*n)] * Y_sim[t - l, :]
        end
        # Constant term
        y_val .+= pop_model.B[:, end]
        y_val .+= shocks[t, :]
        Y_sim[t, :] = y_val
    end
    
    return Y_sim[(burn_in + 1):end, :]
end

end # module
