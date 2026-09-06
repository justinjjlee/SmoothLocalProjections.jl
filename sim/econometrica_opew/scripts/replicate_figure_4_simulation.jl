# scripts/replicate_figure_4_simulation.jl
# Replicates Figure 4 from:
# Montiel Olea, Plagborg-Møller, Qian, and Wolf (Econometrica, 2026)
# "Double Robustness of Local Projections and Some Unpleasant VARithmetic"

using CSV, DataFrames, LinearAlgebra, Statistics, Distributions, Plots, Printf

include(normpath(joinpath(@__DIR__, "..", "src", "Simulation.jl")))
using .Simulation

# Ensure output directory exists
fig_dir = normpath(joinpath(@__DIR__, "..", "figures"))
mkpath(fig_dir)

println("===================================================================")
println("REPLICATING FIGURE 4: EMPIRICALLY CALIBRATED SIMULATION STUDY")
println("Calibration: Känzig (2021, AER) Oil Supply News Application")
println("===================================================================\n")

# 1. Load cleaned empirical dataset
data_csv = normpath(joinpath(@__DIR__, "..", "data", "kanzig_cleaned_data.csv"))
if !isfile(data_csv)
    println("Cleaned dataset not found. Running prepare_data_and_dgp.jl...")
    include("prepare_data_and_dgp.jl")
end

df_clean = CSV.read(data_csv, DataFrame)
data_oil = Matrix(df_clean)
T_raw, n_y = size(data_oil)
println("Loaded empirical dataset: ", size(data_oil), " (7 variables, 528 months)")

# 2. Estimate Population VAR(18) (VAR_pop)
p_pop = 18
T_eff18 = T_raw - p_pop
Y18 = data_oil[(p_pop + 1):end, :]
X18 = zeros(T_eff18, n_y * p_pop)
for l in 1:p_pop
    X18[:, ((l - 1)*n_y + 1):(l*n_y)] = data_oil[(p_pop - l + 1):(end - l), :]
end
B_pop = Matrix((X18 \ Y18)')
resid_pop = Y18 - X18 * B_pop'
Sigma_pop = (resid_pop' * resid_pop) ./ (T_eff18 - size(B_pop, 2))

# 3. Estimate Misspecified Baseline VAR(12) (VAR_p)
p_12 = 12
T_eff12 = T_raw - p_12
Y12 = data_oil[(p_12 + 1):end, :]
X12 = zeros(T_eff12, n_y * p_12)
for l in 1:p_12
    X12[:, ((l - 1)*n_y + 1):(l*n_y)] = data_oil[(p_12 - l + 1):(end - l), :]
end
B_12 = Matrix((X12 \ Y12)')
resid_12 = Y12 - X12 * B_12'
Sigma_12 = (resid_12' * resid_12) ./ (T_eff12 - size(B_12, 2))

# 4. Compute Wold IRFs for VAR_pop and residual VMA
VMA_hor = 200
C_pop = cholesky(Hermitian(Sigma_pop)).L

IRF_pop = zeros(VMA_hor, n_y, n_y)
IRF_pop[1, :, :] = C_pop

for l in 2:VMA_hor
    for j in 1:min(l - 1, p_pop)
        Aj = B_pop[:, ((j - 1)*n_y + 1):(j*n_y)]
        IRF_pop[l, :, :] += Aj * IRF_pop[l - j, :, :]
    end
end

VMA_IRF = zeros(VMA_hor, n_y, n_y)
for l in 1:VMA_hor
    VMA_IRF[l, :, :] += IRF_pop[l, :, :]
    for q in 1:min(l - 1, p_12)
        Aq = -B_12[:, ((q - 1)*n_y + 1):(q*n_y)]
        VMA_IRF[l, :, :] += Aq * IRF_pop[l - q, :, :]
    end
end

# 5. Extract H and D from VMA(1)
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
    A_c[(n_y + 1):end, 1:(n_y*(p_12 - 1))] = Matrix(1.0I, n_y*(p_12 - 1), n_y*(p_12 - 1))
end

H_c = zeros(np, n_y)
H_c[1:n_y, :] = H

# 6. Build Local-to-SVAR ABCD State-Space System
alpha_lags = 100
T_sim = 720
n_s = np + alpha_lags * n_y

A_state = zeros(n_s, n_s)
A_state[1:np, 1:np] = A_c

sqrt_D = Diagonal(sqrt.(D_diag))
for i in 1:alpha_lags
    col_idx = np + (i - 1)*n_y + 1 : np + i*n_y
    A_state[1:np, col_idx] = H_c * (T_sim^(-0.5)) * alpha_tilde[i + 1, :, :] * sqrt_D
end

if alpha_lags > 1
    A_state[(np + n_y + 1):n_s, (np + 1):(n_s - n_y)] = Matrix(1.0I, (alpha_lags - 1)*n_y, (alpha_lags - 1)*n_y)
end

B_state = zeros(n_s, n_y)
B_state[1:np, :] = H_c * sqrt_D
B_state[(np + 1):(np + n_y), :] = Matrix(1.0I, n_y, n_y)

C_obs = A_state[1:n_y, :]
D_obs = B_state[1:n_y, :]

# 7. Compute True Population IRFs of CPI (variable 7) to Oil News Shock (variable 1)
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

println("True population IRFs computed:")
println("   - Horizon 0:  ", round(true_irfs[1], digits=4))
println("   - Horizon 12: ", round(true_irfs[13], digits=4))
println("   - Horizon 24: ", round(true_irfs[25], digits=4))
println("   - Horizon 50: ", round(true_irfs[51], digits=4))

# 8. Data simulation routine from ABCD system
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

# 9. Monte Carlo Simulation Loop
N_reps = 500 # High accuracy replication (500 reps)
z_crit = quantile(Normal(0, 1), 0.95) # 1.64485

println("\nStarting Monte Carlo simulation with N = $N_reps replications (T = $T_sim)...")

cov_var_12 = zeros(H_max + 1)
cov_lp_12  = zeros(H_max + 1)
len_var_12 = [Float64[] for _ in 0:H_max]
len_lp_12  = [Float64[] for _ in 0:H_max]

cov_var_aic = zeros(H_max + 1)
cov_lp_aic  = zeros(H_max + 1)
len_var_aic = [Float64[] for _ in 0:H_max]
len_lp_aic  = [Float64[] for _ in 0:H_max]
aic_lags_recorded = Float64[]

for rep in 1:N_reps
    if rep % 100 == 0 || rep == 1
        @printf("  Progress: %4d / %4d replications\n", rep, N_reps)
    end
    
    Y_sim = sim_abcd(T_sim)
    
    # --- Experiment 1: Fixed p = 12 ---
    v_12 = estimate_var(Y_sim, 12)
    irf_v12, se_v12 = Simulation.compute_var_irf_and_se(v_12, H_max, 7)
    irf_l12, se_l12 = estimate_lp(Y_sim, 12, H_max, 7, 1)
    
    for h in 0:H_max
        # VAR(12)
        v_low = irf_v12[h + 1] - z_crit * se_v12[h + 1]
        v_upp = irf_v12[h + 1] + z_crit * se_v12[h + 1]
        if v_low <= true_irfs[h + 1] <= v_upp
            cov_var_12[h + 1] += 1.0
        end
        push!(len_var_12[h + 1], v_upp - v_low)
        
        # LP(12)
        l_low = irf_l12[h + 1] - z_crit * se_l12[h + 1]
        l_upp = irf_l12[h + 1] + z_crit * se_l12[h + 1]
        if l_low <= true_irfs[h + 1] <= l_upp
            cov_lp_12[h + 1] += 1.0
        end
        push!(len_lp_12[h + 1], l_upp - l_low)
    end
    
    # --- Experiment 2: Lag selection via AIC (p_max = 24) ---
    p_aic = select_lag_aic(Y_sim, 24)
    push!(aic_lags_recorded, Float64(p_aic))
    
    v_aic = estimate_var(Y_sim, p_aic)
    irf_vaic, se_vaic = Simulation.compute_var_irf_and_se(v_aic, H_max, 7)
    irf_laic, se_laic = estimate_lp(Y_sim, p_aic, H_max, 7, 1)
    
    for h in 0:H_max
        # VAR(AIC)
        va_low = irf_vaic[h + 1] - z_crit * se_vaic[h + 1]
        va_upp = irf_vaic[h + 1] + z_crit * se_vaic[h + 1]
        if va_low <= true_irfs[h + 1] <= va_upp
            cov_var_aic[h + 1] += 1.0
        end
        push!(len_var_aic[h + 1], va_upp - va_low)
        
        # LP(AIC)
        la_low = irf_laic[h + 1] - z_crit * se_laic[h + 1]
        la_upp = irf_laic[h + 1] + z_crit * se_laic[h + 1]
        if la_low <= true_irfs[h + 1] <= la_upp
            cov_lp_aic[h + 1] += 1.0
        end
        push!(len_lp_aic[h + 1], la_upp - la_low)
    end
end

# Compute empirical rates
cov_var_12 ./= N_reps
cov_lp_12  ./= N_reps
cov_var_aic ./= N_reps
cov_lp_aic  ./= N_reps

med_len_var_12 = [median(len_var_12[h + 1]) for h in 0:H_max]
med_len_lp_12  = [median(len_lp_12[h + 1]) for h in 0:H_max]

med_len_var_aic = [median(len_var_aic[h + 1]) for h in 0:H_max]
med_len_lp_aic  = [median(len_lp_aic[h + 1]) for h in 0:H_max]

mean_aic_lag = mean(aic_lags_recorded)

println("\n===================================================================")
println("SIMULATION RESULTS BENCHMARKS VS. PUBLISHED ECONOMETRICA PAPER")
println("===================================================================")
@printf("1. Mean lag length selected by AIC: %.2f (Paper text: 9.7)\n", mean_aic_lag)
@printf("2. Horizon 0 coverage: VAR(12) = %.1f%% | LP(12) = %.1f%%\n", cov_var_12[1]*100, cov_lp_12[1]*100)
@printf("3. Horizon 25 coverage (p = 12):\n")
@printf("   - VAR(12): %.1f%% (Paper text: falls below 60%%)\n", cov_var_12[26]*100)
@printf("   - LP(12):  %.1f%% (Paper text: ~90%% nominal)\n", cov_lp_12[26]*100)
@printf("4. Horizon 50 coverage (p = 12):\n")
@printf("   - VAR(12): %.1f%% (Paper text: ~50%%)\n", cov_var_12[51]*100)
@printf("   - LP(12):  %.1f%% (Paper text: ~90%%)\n", cov_lp_12[51]*100)
@printf("5. Horizon 50 coverage (AIC):\n")
@printf("   - VAR(AIC): %.1f%% (Paper text: falls below 60%%)\n", cov_var_aic[51]*100)
@printf("   - LP(AIC):  %.1f%% (Paper text: ~90%%)\n", cov_lp_aic[51]*100)
println("===================================================================")

# -----------------------------------------------------------------------------
# PLOT FIGURE 4 (4 SUBPLOTS MATCHING ECONOMETRICA 2026 PAGE 20)
# -----------------------------------------------------------------------------
horizons = 0:H_max

Plots.default(
    fontfamily = "sans-serif",
    titlefontsize = 10,
    guidefontsize = 9,
    tickfontsize = 8,
    legendfontsize = 8,
    linewidth = 2.0,
    dpi = 300
)

# Top Left: Coverage for p = 12
p_tl = plot(
    horizons, cov_var_12,
    color = :red, linewidth = 2.2, label = "VAR",
    title = "LAG LENGTH p = 12\ncoverage probability",
    xlabel = "horizon", ylabel = "",
    xlims = (0, 50), ylims = (0.0, 1.0),
    legend = :bottomright
)
plot!(p_tl, horizons, cov_lp_12, color = :blue, linestyle = :dot, linewidth = 2.2, label = "LP")
hline!(p_tl, [0.90], color = :black, linestyle = :solid, linewidth = 0.8, label = false)

# Top Right: Median length, log scale for p = 12
p_tr = plot(
    horizons, med_len_var_12,
    color = :red, linewidth = 2.2, label = "VAR",
    title = "LAG LENGTH p = 12\nmedian length, log scale",
    xlabel = "horizon", ylabel = "",
    yscale = :log10,
    xlims = (0, 50), ylims = (1e-3, 1e1),
    yticks = ([1e-3, 1e-2, 1e-1, 1e0, 1e1], ["0.001", "0.01", "0.1", "1", "10"]),
    legend = false
)
plot!(p_tr, horizons, med_len_lp_12, color = :blue, linestyle = :dot, linewidth = 2.2, label = "LP")

# Bottom Left: Coverage for AIC
p_bl = plot(
    horizons, cov_var_aic,
    color = :red, linewidth = 2.2, label = "VAR",
    title = "LAG LENGTH VIA AIC\ncoverage probability",
    xlabel = "horizon", ylabel = "",
    xlims = (0, 50), ylims = (0.0, 1.0),
    legend = :bottomright
)
plot!(p_bl, horizons, cov_lp_aic, color = :blue, linestyle = :dot, linewidth = 2.2, label = "LP")
hline!(p_bl, [0.90], color = :black, linestyle = :solid, linewidth = 0.8, label = false)

# Bottom Right: Median length, log scale for AIC
p_br = plot(
    horizons, med_len_var_aic,
    color = :red, linewidth = 2.2, label = "VAR",
    title = "LAG LENGTH VIA AIC\nmedian length, log scale",
    xlabel = "horizon", ylabel = "",
    yscale = :log10,
    xlims = (0, 50), ylims = (1e-3, 1e1),
    yticks = ([1e-3, 1e-2, 1e-1, 1e0, 1e1], ["0.001", "0.01", "0.1", "1", "10"]),
    legend = false
)
plot!(p_br, horizons, med_len_lp_aic, color = :blue, linestyle = :dot, linewidth = 2.2, label = "LP")

p_fig4 = plot(p_tl, p_tr, p_bl, p_br, layout = (2, 2), size = (900, 750))
fig4_path = joinpath(fig_dir, "figure_4_simulation.png")
savefig(p_fig4, fig4_path)
println("Successfully saved Figure 4 to: ", fig4_path)
