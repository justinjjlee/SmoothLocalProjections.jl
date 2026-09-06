# scripts/test_varma_simulation.jl
using CSV, DataFrames, LinearAlgebra, Statistics, Distributions, Printf

df = CSV.read(normpath(joinpath(@__DIR__, "..", "data", "kanzig_cleaned_data.csv")), DataFrame)
data_oil = Matrix(df)
T_raw, n_y = size(data_oil)

# 1. Estimate VAR(18) (VAR_pop)
p_pop = 18
T_eff18 = T_raw - p_pop
Y18 = data_oil[(p_pop+1):end, :]
X18 = zeros(T_eff18, n_y * p_pop)
for l in 1:p_pop
    X18[:, ((l-1)*n_y + 1):(l*n_y)] = data_oil[(p_pop - l + 1):(end - l), :]
end
B_pop = Matrix((X18 \ Y18)')
resid_pop = Y18 - X18 * B_pop'
Sigma_pop = (resid_pop' * resid_pop) ./ (T_eff18 - size(B_pop, 2))

# 2. Estimate VAR(12) (VAR_p)
p_12 = 12
T_eff12 = T_raw - p_12
Y12 = data_oil[(p_12+1):end, :]
X12 = zeros(T_eff12, n_y * p_12)
for l in 1:p_12
    X12[:, ((l-1)*n_y + 1):(l*n_y)] = data_oil[(p_12 - l + 1):(end - l), :]
end
B_12 = Matrix((X12 \ Y12)')
resid_12 = Y12 - X12 * B_12'
Sigma_12 = (resid_12' * resid_12) ./ (T_eff12 - size(B_12, 2))

# 3. Compute Wold IRFs for VAR_pop and residual VMA
VMA_hor = 200
C_pop = cholesky(Hermitian(Sigma_pop)).L

IRF_pop = zeros(VMA_hor, n_y, n_y)
IRF_pop[1, :, :] = C_pop

for l in 2:VMA_hor
    for j in 1:min(l - 1, p_pop)
        Aj = B_pop[:, ((j-1)*n_y + 1):(j*n_y)]
        IRF_pop[l, :, :] += Aj * IRF_pop[l - j, :, :]
    end
end

VMA_IRF = zeros(VMA_hor, n_y, n_y)
for l in 1:VMA_hor
    VMA_IRF[l, :, :] += IRF_pop[l, :, :]
    for q in 1:min(l - 1, p_12)
        Aq = -B_12[:, ((q-1)*n_y + 1):(q*n_y)]
        VMA_IRF[l, :, :] += Aq * IRF_pop[l - q, :, :]
    end
end

# 4. H and D from VMA(1)
HD = cholesky(Hermitian(VMA_IRF[1, :, :] * VMA_IRF[1, :, :]')).L
H = zeros(n_y, n_y)
D_diag = zeros(n_y)
for i in 1:n_y
    H[:, i] = HD[:, i] ./ HD[i, i]
    D_diag[i] = HD[i, i]^2
end

alpha_tilde = zeros(VMA_hor, n_y, n_y)
inv_H = inv(H)
inv_VMA1 = inv(VMA_IRF[1, :, :])
for l in 1:VMA_hor
    alpha_tilde[l, :, :] = inv_H * VMA_IRF[l, :, :] * inv_VMA1 * H
end

# Companion matrices
np = n_y * p_12
A_c = zeros(np, np)
A_c[1:n_y, :] = B_12
if p_12 > 1
    A_c[(n_y+1):end, 1:(n_y*(p_12-1))] = Matrix(1.0I, n_y*(p_12-1), n_y*(p_12-1))
end

H_c = zeros(np, n_y)
H_c[1:n_y, :] = H

# 5. Build State Space ABCD system
alpha_lags = 100
T_sim = 720
zeta = 0.5
n_s = np + alpha_lags * n_y

A_state = zeros(n_s, n_s)
A_state[1:np, 1:np] = A_c

sqrt_D = Diagonal(sqrt.(D_diag))
for i in 1:alpha_lags
    # In simul_oil.m: alpha_tilde(2:end) was scaled by T_scale^zeta, then divided by T^zeta in set_up_varma
    # So net multiplier is 1.0!
    col_idx = np + (i - 1)*n_y + 1 : np + i*n_y
    A_state[1:np, col_idx] = H_c * alpha_tilde[i + 1, :, :] * sqrt_D
end

if alpha_lags > 1
    A_state[(np + n_y + 1):n_s, (np + 1):(n_s - n_y)] = Matrix(1.0I, (alpha_lags - 1)*n_y, (alpha_lags - 1)*n_y)
end

B_state = zeros(n_s, n_y)
B_state[1:np, :] = H_c * sqrt_D
B_state[(np + 1):(np + n_y), :] = Matrix(1.0I, n_y, n_y)

C_obs = A_state[1:n_y, :]
D_obs = B_state[1:n_y, :]

# 6. Compute true population IRFs via ABCD system
H_max = 50
shock_weight = zeros(n_y)
shock_weight[1] = 1.0 / (H[1, 1] * sqrt(D_diag[1]))

true_irfs = zeros(H_max + 1)
# h = 0
true_irfs[1] = dot(D_obs[7, :], shock_weight)

let
    global A_pow = Matrix(1.0I, n_s, n_s)
    for h in 1:H_max
        resp = C_obs * A_pow * B_state * shock_weight
        true_irfs[h + 1] = resp[7]
        A_pow = A_pow * A_state
    end
end

println("True IRFs from VARMA ABCD system:")
println("  h = 0:  ", round(true_irfs[1], digits=4))
println("  h = 12: ", round(true_irfs[13], digits=4))
println("  h = 24: ", round(true_irfs[25], digits=4))
println("  h = 50: ", round(true_irfs[51], digits=4))
