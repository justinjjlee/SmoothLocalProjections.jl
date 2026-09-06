# scripts/test_simulation_quick.jl
using CSV, DataFrames, LinearAlgebra, Statistics, Distributions, Printf

include(normpath(joinpath(@__DIR__, "..", "src", "Simulation.jl")))
using .Simulation

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
n_s = np + alpha_lags * n_y

A_state = zeros(n_s, n_s)
A_state[1:np, 1:np] = A_c

sqrt_D = Diagonal(sqrt.(D_diag))
for i in 1:alpha_lags
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

# 6. Compute true population IRFs
H_max = 50
shock_weight = zeros(n_y)
shock_weight[1] = 1.0 / (H[1, 1] * sqrt(D_diag[1]))

true_irfs = zeros(H_max + 1)
true_irfs[1] = dot(D_obs[7, :], shock_weight)
let
    global A_pow = Matrix(1.0I, n_s, n_s)
    for h in 1:H_max
        resp = C_obs * A_pow * B_state * shock_weight
        true_irfs[h + 1] = resp[7]
        A_pow = A_pow * A_state
    end
end

# 7. Simulate data function
function sim_abcd(T_len::Int)
    eps = randn(T_len, n_y)
    s = zeros(n_s)
    data_s = zeros(T_len, n_s)
    for t in 1:T_len
        s = A_state * s + B_state * eps[t, :]
        data_s[t, :] = s
    end
    y_out = zeros(T_len, n_y)
    y_out[1, :] = D_obs * eps[1, :]
    for t in 2:T_len
        y_out[t, :] = C_obs * data_s[t - 1, :] + D_obs * eps[t, :]
    end
    return y_out
end

# Run 100 test replications
N_reps = 100
z_crit = 1.64485
cov_var_12 = zeros(H_max + 1)
cov_lp_12 = zeros(H_max + 1)

for rep in 1:N_reps
    Y_sim = sim_abcd(T_sim)
    v_mod = estimate_var(Y_sim, 12)
    irf_v, se_v = Simulation.compute_var_irf_and_se(v_mod, H_max, 7)
    irf_l, se_l = estimate_lp(Y_sim, 12, H_max, 7, 1)
    
    for h in 0:H_max
        if (irf_v[h+1] - z_crit*se_v[h+1]) <= true_irfs[h+1] <= (irf_v[h+1] + z_crit*se_v[h+1])
            cov_var_12[h+1] += 1.0
        end
        if (irf_l[h+1] - z_crit*se_l[h+1]) <= true_irfs[h+1] <= (irf_l[h+1] + z_crit*se_l[h+1])
            cov_lp_12[h+1] += 1.0
        end
    end
end

cov_var_12 ./= N_reps
cov_lp_12 ./= N_reps

println("Coverage check for p = 12:")
println("  h = 0:  VAR = $(round(cov_var_12[1]*100, digits=1))% | LP = $(round(cov_lp_12[1]*100, digits=1))%")
println("  h = 10: VAR = $(round(cov_var_12[11]*100, digits=1))% | LP = $(round(cov_lp_12[11]*100, digits=1))%")
println("  h = 20: VAR = $(round(cov_var_12[21]*100, digits=1))% | LP = $(round(cov_lp_12[21]*100, digits=1))%")
println("  h = 30: VAR = $(round(cov_var_12[31]*100, digits=1))% | LP = $(round(cov_lp_12[31]*100, digits=1))%")
println("  h = 50: VAR = $(round(cov_var_12[51]*100, digits=1))% | LP = $(round(cov_lp_12[51]*100, digits=1))%")
