# scripts/test_delta_se.jl
using LinearAlgebra, CSV, DataFrames

df = CSV.read(normpath(joinpath(@__DIR__, "..", "data", "kanzig_cleaned_data.csv")), DataFrame)
data_oil = Matrix(df)
T_raw, n = size(data_oil)
p = 12
T_eff = T_raw - p
Y_lhs = data_oil[(p+1):end, :]
X = zeros(T_eff, n * p)
for l in 1:p
    X[:, ((l-1)*n + 1):(l*n)] = data_oil[(p - l + 1):(end - l), :]
end
X_expand = hcat(X, ones(T_eff, 1))
betahat = Matrix((X_expand \ Y_lhs)')
B = betahat[:, 1:end-1]
resid = Y_lhs - X_expand * betahat'
Sigma_u = (resid' * resid) ./ (T_eff - size(betahat, 2))

C_chol = cholesky(Hermitian(Sigma_u)).L
nu = C_chol[:, 1] ./ C_chol[1, 1]

np = n * p
A_comp = zeros(np, np)
A_comp[1:n, :] = B
if p > 1
    A_comp[(n+1):end, 1:(n*(p-1))] = Matrix(1.0I, n*(p-1), n*(p-1))
end

nu_comp = zeros(np)
nu_comp[1:n] = nu

e_i = zeros(np); e_i[7] = 1.0

A_pows = [Matrix(1.0I, np, np)]
for h in 1:50
    push!(A_pows, A_pows[end] * A_comp)
end

X_slopes = X_expand[:, 1:end-1]
inv_XX = inv(X_slopes' * X_slopes)

println("Delta method SEs for CPI to oil shock:")
for h in [1, 5, 12, 24, 50]
    Psi_h = zeros(np, np)
    for l in 1:h
        Psi_h .+= A_pows[h - l + 1] * nu_comp * (e_i' * A_pows[l])
    end
    G = Matrix(Psi_h[:, 1:n]') # n × np
    var_h = tr(G * inv_XX * G' * Sigma_u)
    irf_h = dot(e_i, A_pows[h+1] * nu_comp)
    se_h = sqrt(var_h)
    width_h = 2.0 * 1.64485 * se_h
    println("  h = $h: irf = $(round(irf_h, digits=4)), se = $(round(se_h, digits=4)), 90% width = $(round(width_h, digits=4))")
end
