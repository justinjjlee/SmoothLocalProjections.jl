# scripts/plot_simulation_figures.jl
# Replicates Figure 4 (Main Text) and Figure D.1 (Appendix) from:
# Montiel Olea, Plagborg-Møller, Qian, and Wolf (Econometrica, 2026)
# "Double Robustness of Local Projections and Some Unpleasant VARithmetic"

using CSV
using DataFrames
using Plots
using Printf

fig_dir = normpath(joinpath(@__DIR__, "..", "figures"))
data_dir = normpath(joinpath(@__DIR__, "..", "data"))
mkpath(fig_dir)

println("===================================================================")
println("PLOTTING EMPIRICAL SIMULATION RESULTS (KÄNZIG 2021 CALIBRATION)")
println("===================================================================\n")

# Load exported simulation datasets
df_fig4 = CSV.read(joinpath(data_dir, "sim_figure4_data.csv"), DataFrame)
df_d1   = CSV.read(joinpath(data_dir, "sim_figureD1_data.csv"), DataFrame)

# Common styling
Plots.default(
    fontfamily = "sans-serif",
    titlefontsize = 11,
    guidefontsize = 10,
    tickfontsize = 9,
    legendfontsize = 8,
    linewidth = 2.2,
    dpi = 300
)

# Colors and line styles matching Econometrica paper / plot_oil.m
c_var    = RGB(204/255, 0/255, 0/255)         # Strong red
c_var_b  = RGB(230/255, 128/255, 128/255)     # Light red
c_lp     = RGB(102/255, 178/255, 255/255)     # Light blue
c_lp_b   = RGB(51/255, 102/255, 204/255)      # Deep blue

# -----------------------------------------------------------------------------
# FIGURE 4: Coverage and Median Length for p = 12 (Top) and AIC (Bottom)
# -----------------------------------------------------------------------------
df_p12 = filter(row -> row.lag_spec == "p=12", df_fig4)
df_aic = filter(row -> row.lag_spec == "AIC", df_fig4)

horizons = df_p12.horizon

# Panel (1, 1): p=12 Coverage
p11 = plot(
    horizons, df_p12.var_cov,
    color = c_var, linestyle = :solid, label = "VAR",
    title = "LAG LENGTH p = 12\ncoverage probability",
    xlabel = "horizon", ylabel = "",
    xlims = (0, 50), ylims = (0.0, 1.0),
    legend = :bottomright
)
plot!(p11, horizons, df_p12.var_b_cov, color = c_var_b, linestyle = :dashdot, label = "VAR_b")
plot!(p11, horizons, df_p12.lp_cov,    color = c_lp,    linestyle = :dot,     label = "LP")
plot!(p11, horizons, df_p12.lp_b_cov,  color = c_lp_b,  linestyle = :dash,    label = "LP_b")
hline!(p11, [0.90], color = :black, linestyle = :dot, linewidth = 1.0, label = false)

# Panel (1, 2): p=12 Median Length (log scale)
p12 = plot(
    horizons, df_p12.var_len,
    color = c_var, linestyle = :solid, label = "VAR",
    title = "LAG LENGTH p = 12\nmedian length, log scale",
    xlabel = "horizon", ylabel = "",
    yscale = :log10,
    xlims = (0, 50), ylims = (1e-3, 1e1),
    yticks = ([1e-3, 1e-2, 1e-1, 1e0, 1e1], ["0.001", "0.01", "0.1", "1", "10"]),
    legend = false
)
plot!(p12, horizons, df_p12.var_b_len, color = c_var_b, linestyle = :dashdot, label = "VAR_b")
plot!(p12, horizons, df_p12.lp_len,    color = c_lp,    linestyle = :dot,     label = "LP")
plot!(p12, horizons, df_p12.lp_b_len,  color = c_lp_b,  linestyle = :dash,    label = "LP_b")

# Panel (2, 1): AIC Coverage
p21 = plot(
    horizons, df_aic.var_cov,
    color = c_var, linestyle = :solid, label = "VAR",
    title = "LAG LENGTH VIA AIC\ncoverage probability",
    xlabel = "horizon", ylabel = "",
    xlims = (0, 50), ylims = (0.0, 1.0),
    legend = :bottomright
)
plot!(p21, horizons, df_aic.var_b_cov, color = c_var_b, linestyle = :dashdot, label = "VAR_b")
plot!(p21, horizons, df_aic.lp_cov,    color = c_lp,    linestyle = :dot,     label = "LP")
plot!(p21, horizons, df_aic.lp_b_cov,  color = c_lp_b,  linestyle = :dash,    label = "LP_b")
hline!(p21, [0.90], color = :black, linestyle = :dot, linewidth = 1.0, label = false)

# Panel (2, 2): AIC Median Length (log scale)
p22 = plot(
    horizons, df_aic.var_len,
    color = c_var, linestyle = :solid, label = "VAR",
    title = "LAG LENGTH VIA AIC\nmedian length, log scale",
    xlabel = "horizon", ylabel = "",
    yscale = :log10,
    xlims = (0, 50), ylims = (1e-3, 1e1),
    yticks = ([1e-3, 1e-2, 1e-1, 1e0, 1e1], ["0.001", "0.01", "0.1", "1", "10"]),
    legend = false
)
plot!(p22, horizons, df_aic.var_b_len, color = c_var_b, linestyle = :dashdot, label = "VAR_b")
plot!(p22, horizons, df_aic.lp_len,    color = c_lp,    linestyle = :dot,     label = "LP")
plot!(p22, horizons, df_aic.lp_b_len,  color = c_lp_b,  linestyle = :dash,    label = "LP_b")

p_fig4 = plot(p11, p12, p21, p22, layout = (2, 2), size = (920, 760))
fig4_out = joinpath(fig_dir, "figure_4_simulation.png")
savefig(p_fig4, fig4_out)
println("[1/2] Successfully generated Figure 4: ", fig4_out)

# -----------------------------------------------------------------------------
# FIGURE D.1: Simulation Results for p = 15 and p = 18 (Online Appendix)
# -----------------------------------------------------------------------------
df_p15 = filter(row -> row.lag_spec == "p=15", df_d1)
df_p18 = filter(row -> row.lag_spec == "p=18", df_d1)

# p=15 Coverage & Length
p_d1_11 = plot(
    horizons, df_p15.var_cov, color = c_var, linestyle = :solid, label = "VAR",
    title = "LAG LENGTH p = 15\ncoverage probability", xlabel = "horizon", ylabel = "",
    xlims = (0, 50), ylims = (0.0, 1.0), legend = :bottomright
)
plot!(p_d1_11, horizons, df_p15.var_b_cov, color = c_var_b, linestyle = :dashdot, label = "VAR_b")
plot!(p_d1_11, horizons, df_p15.lp_cov,    color = c_lp,    linestyle = :dot,     label = "LP")
plot!(p_d1_11, horizons, df_p15.lp_b_cov,  color = c_lp_b,  linestyle = :dash,    label = "LP_b")
hline!(p_d1_11, [0.90], color = :black, linestyle = :dot, linewidth = 1.0, label = false)

p_d1_12 = plot(
    horizons, df_p15.var_len, color = c_var, linestyle = :solid, label = "VAR",
    title = "LAG LENGTH p = 15\nmedian length, log scale", xlabel = "horizon", ylabel = "",
    yscale = :log10, xlims = (0, 50), ylims = (1e-3, 1e1),
    yticks = ([1e-3, 1e-2, 1e-1, 1e0, 1e1], ["0.001", "0.01", "0.1", "1", "10"]), legend = false
)
plot!(p_d1_12, horizons, df_p15.var_b_len, color = c_var_b, linestyle = :dashdot, label = "VAR_b")
plot!(p_d1_12, horizons, df_p15.lp_len,    color = c_lp,    linestyle = :dot,     label = "LP")
plot!(p_d1_12, horizons, df_p15.lp_b_len,  color = c_lp_b,  linestyle = :dash,    label = "LP_b")

# p=18 Coverage & Length
p_d1_21 = plot(
    horizons, df_p18.var_cov, color = c_var, linestyle = :solid, label = "VAR",
    title = "LAG LENGTH p = 18\ncoverage probability", xlabel = "horizon", ylabel = "",
    xlims = (0, 50), ylims = (0.0, 1.0), legend = :bottomright
)
plot!(p_d1_21, horizons, df_p18.var_b_cov, color = c_var_b, linestyle = :dashdot, label = "VAR_b")
plot!(p_d1_21, horizons, df_p18.lp_cov,    color = c_lp,    linestyle = :dot,     label = "LP")
plot!(p_d1_21, horizons, df_p18.lp_b_cov,  color = c_lp_b,  linestyle = :dash,    label = "LP_b")
hline!(p_d1_21, [0.90], color = :black, linestyle = :dot, linewidth = 1.0, label = false)

p_d1_22 = plot(
    horizons, df_p18.var_len, color = c_var, linestyle = :solid, label = "VAR",
    title = "LAG LENGTH p = 18\nmedian length, log scale", xlabel = "horizon", ylabel = "",
    yscale = :log10, xlims = (0, 50), ylims = (1e-3, 1e1),
    yticks = ([1e-3, 1e-2, 1e-1, 1e0, 1e1], ["0.001", "0.01", "0.1", "1", "10"]), legend = false
)
plot!(p_d1_22, horizons, df_p18.var_b_len, color = c_var_b, linestyle = :dashdot, label = "VAR_b")
plot!(p_d1_22, horizons, df_p18.lp_len,    color = c_lp,    linestyle = :dot,     label = "LP")
plot!(p_d1_22, horizons, df_p18.lp_b_len,  color = c_lp_b,  linestyle = :dash,    label = "LP_b")

p_figD1 = plot(p_d1_11, p_d1_12, p_d1_21, p_d1_22, layout = (2, 2), size = (920, 760))
figD1_out = joinpath(fig_dir, "figure_d1_simulation_appendix.png")
savefig(p_figD1, figD1_out)
println("[2/2] Successfully generated Figure D.1: ", figD1_out)

println("\n===================================================================")
println("VALIDATION OF PUBLISHED NUMERICAL METRICS")
println("===================================================================")
println("1. Fixed Lag p = 12 Coverage:")
@printf("   - h = 0:  VAR = %.1f%%, VAR_b = %.1f%%, LP = %.1f%%, LP_b = %.1f%%\n",
        df_p12.var_cov[1]*100, df_p12.var_b_cov[1]*100, df_p12.lp_cov[1]*100, df_p12.lp_b_cov[1]*100)
@printf("   - h = 25: VAR = %.1f%%, VAR_b = %.1f%%, LP = %.1f%%, LP_b = %.1f%%\n",
        df_p12.var_cov[26]*100, df_p12.var_b_cov[26]*100, df_p12.lp_cov[26]*100, df_p12.lp_b_cov[26]*100)
@printf("   - h = 50: VAR = %.1f%%, VAR_b = %.1f%%, LP = %.1f%%, LP_b = %.1f%%\n",
        df_p12.var_cov[51]*100, df_p12.var_b_cov[51]*100, df_p12.lp_cov[51]*100, df_p12.lp_b_cov[51]*100)

println("2. AIC Lag Selection Coverage:")
@printf("   - h = 0:  VAR = %.1f%%, VAR_b = %.1f%%, LP = %.1f%%, LP_b = %.1f%%\n",
        df_aic.var_cov[1]*100, df_aic.var_b_cov[1]*100, df_aic.lp_cov[1]*100, df_aic.lp_b_cov[1]*100)
@printf("   - h = 40: VAR = %.1f%%, VAR_b = %.1f%%, LP = %.1f%%, LP_b = %.1f%%\n",
        df_aic.var_cov[41]*100, df_aic.var_b_cov[41]*100, df_aic.lp_cov[41]*100, df_aic.lp_b_cov[41]*100)
@printf("   - h = 50: VAR = %.1f%%, VAR_b = %.1f%%, LP = %.1f%%, LP_b = %.1f%%\n",
        df_aic.var_cov[51]*100, df_aic.var_b_cov[51]*100, df_aic.lp_cov[51]*100, df_aic.lp_b_cov[51]*100)

println("3. Median Length (p = 12):")
@printf("   - h = 0:  VAR = %.4f, LP = %.4f (ratio: %.3f)\n",
        df_p12.var_len[1], df_p12.lp_len[1], df_p12.var_len[1]/df_p12.lp_len[1])
@printf("   - h = 50: VAR = %.4f, LP = %.4f (ratio: %.3f)\n",
        df_p12.var_len[51], df_p12.lp_len[51], df_p12.var_len[51]/df_p12.lp_len[51])
println("===================================================================")
