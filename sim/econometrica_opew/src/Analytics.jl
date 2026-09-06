module Analytics

using Distributions
using LinearAlgebra
using Statistics

export r_tail, solve_cv, worst_case_coverage, worst_case_joint_prob
export bias_aware_rel_length, least_favorable_alpha_ar1, optimal_averaging_ci

"""
    r_tail(b::Float64, c::Float64)

Compute two-tailed folded normal probability r(b; c) = P(|Z + b| > c),
where Z ~ N(0, 1). Defined in Corollary 3.1:
r(b; c) = Φ(-c - b) + Φ(-c + b).
"""
function r_tail(b::Float64, c::Float64)
    d = Normal(0, 1)
    return cdf(d, -c - b) + cdf(d, -c + b)
end

"""
    solve_cv(b::Float64, a::Float64 = 0.10)

Find the bias-aware critical value cv_{1-a}(b) solving:
r(b; cv) = a.
Uses robust bisection search.
"""
function solve_cv(b::Float64, a::Float64 = 0.10)
    z_init = quantile(Normal(0, 1), 1.0 - a / 2.0)
    if abs(b) < 1e-10
        return z_init
    end
    low = z_init
    high = z_init + abs(b) + 5.0
    for _ in 1:60
        mid = 0.5 * (low + high)
        if r_tail(b, mid) > a
            low = mid
        else
            high = mid
        end
    end
    return 0.5 * (low + high)
end

"""
    worst_case_coverage(R::Float64, M::Float64, a::Float64 = 0.10)

Evaluate the worst-case asymptotic coverage probability of the conventional
level-(1 - a) SVAR confidence interval (Corollary 4.3):
inf P(θ ∈ CI(δ)) = 1 - r(M * τ; z_{1-a/2}),
where τ = √(1 / R^2 - 1) and R = √(aVar(δ) / aVar(β)).
"""
function worst_case_coverage(R::Float64, M::Float64, a::Float64 = 0.10)
    z = quantile(Normal(0, 1), 1.0 - a / 2.0)
    if R >= 0.9999
        return 1.0 - a
    elseif R <= 1e-6
        return 0.0
    end
    tau = sqrt(1.0 / (R^2) - 1.0)
    b_star = M * tau
    return 1.0 - r_tail(b_star, z)
end

"""
    worst_case_joint_prob(R::Float64, a::Float64 = 0.10)

Evaluate the worst-case asymptotic probability of the joint event that the
conventional VAR CI fails to cover the true impulse response AND the Hausman
test fails to reject misspecification (Corollary 4.4):
sup_{b ≥ 0} r(b; z_{1-a/2}) * [1 - r(b / τ; z_{1-a/2})].
"""
function worst_case_joint_prob(R::Float64, a::Float64 = 0.10)
    z = quantile(Normal(0, 1), 1.0 - a / 2.0)
    if R >= 0.9999
        return a * (1.0 - a)
    end
    tau = sqrt(1.0 / (R^2) - 1.0)
    
    obj(b) = r_tail(b, z) * (1.0 - r_tail(b / tau, z))
    
    # 2-stage grid search and refinement over b >= 0
    best_val = a * (1.0 - a)
    best_b = 0.0
    for b in 0.0:0.02:8.0
        val = obj(b)
        if val > best_val
            best_val = val
            best_b = b
        end
    end
    # Local refinement
    low_b = max(0.0, best_b - 0.05)
    high_b = best_b + 0.05
    for b in range(low_b, high_b, length=100)
        val = obj(b)
        if val > best_val
            best_val = val
        end
    end
    return best_val
end

"""
    bias_aware_rel_length(R::Float64, M::Float64, a::Float64 = 0.10)

Compute the relative length of the bias-aware VAR confidence interval vs. the
conventional LP interval (Section 4.3):
[cv_{1-a}(M * τ) / z_{1-a/2}] * R.
"""
function bias_aware_rel_length(R::Float64, M::Float64, a::Float64 = 0.10)
    z = quantile(Normal(0, 1), 1.0 - a / 2.0)
    if R >= 0.9999
        return 1.0
    elseif R <= 1e-6
        # As R -> 0, cv(M/R) * R -> M, so rel_len -> M / z
        return M / z
    end
    tau = sqrt(1.0 / (R^2) - 1.0)
    b_star = M * tau
    cv = solve_cv(b_star, a)
    return (cv / z) * R
end

"""
    least_favorable_alpha_ar1(rho::Float64, h::Int, max_lag::Int = 20)

Compute the least-favorable MA polynomial coefficients α†_{ℓ, h}
for a univariate local-to-AR(1) model (Equation (4.1) / Appendix A.1 / plot_arbias.m):
α†_ℓ = h * ρ^{h-1} * (1 - ρ^2) * ρ^{ℓ-1} - 1(ℓ ≤ h) * ρ^{h-ℓ}.
"""
function least_favorable_alpha_ar1(rho::Float64, h::Int, max_lag::Int = 20)
    alpha = zeros(max_lag)
    for l in 1:max_lag
        term1 = h * (rho^(h - 1)) * (1.0 - rho^2) * (rho^(l - 1))
        term2 = (l <= h) ? rho^(h - l) : 0.0
        alpha[l] = term1 - term2
    end
    return alpha
end

"""
    optimal_averaging_ci(R::Float64, M::Float64, a::Float64 = 0.10)

Find the length-optimal weight ω* on LP in the bias-aware model averaging
confidence interval centered at θ̂_h(ω) = ω β̂_h + (1 - ω) δ̂_h (Appendix A.2):
ω* = argmin_{ω ∈ [0, 1]} cv_{1-a}( (1-ω)Mτ / √(1 + ω²τ²) ) * √(1 + ω²τ²).
Returns (omega_star, rel_length).
"""
function optimal_averaging_ci(R::Float64, M::Float64, a::Float64 = 0.10)
    z = quantile(Normal(0, 1), 1.0 - a / 2.0)
    if R >= 0.9999
        return 1.0, 1.0
    end
    tau = sqrt(1.0 / (R^2) - 1.0)
    
    obj(w) = solve_cv((1.0 - w) * M * tau / sqrt(1.0 + w^2 * tau^2), a) * sqrt(1.0 + w^2 * tau^2)
    
    # 1D search over omega in [0, 1]
    best_w = 1.0
    best_len = obj(1.0)
    for w in 0.0:0.01:1.0
        val = obj(w)
        if val < best_len
            best_len = val
            best_w = w
        end
    end
    rel_len = (best_len / z) * R
    return best_w, rel_len
end

end # module
