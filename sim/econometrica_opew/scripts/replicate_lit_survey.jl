# scripts/replicate_lit_survey.jl
# Replicates the empirical literature review in Section 5.1 of:
# Montiel Olea, Plagborg-Møller, Qian, and Wolf (Econometrica, 2026)

using CSV, DataFrames, Statistics, Plots, Printf

data_path = normpath(joinpath(@__DIR__, "..", "data", "lit_varlags_raw.csv"))
df_lit = CSV.read(data_path, DataFrame)

lags = Float64.(df_lit.lags)
maxhor = Float64.(df_lit.maxhor)
freq = Float64.(df_lit.freq)
lag_crit = Float64.(df_lit.lag_crit)
bayes = Float64.(df_lit.bayes)

n_papers = nrow(df_lit)
lags_freq_rel = lags ./ freq
lags_maxhor_rel = lags ./ maxhor

# Statistics cited in Section 5.1 of Econometrica (2026):
avg_lags_freq = mean(lags_freq_rel)
ic_lags_freq = sum(lags_freq_rel .* lag_crit) / sum(lag_crit)

avg_lags_maxhor = mean(lags_maxhor_rel)
ic_lags_maxhor = sum(lags_maxhor_rel .* lag_crit) / sum(lag_crit)

pct_ic = mean(lag_crit) * 100.0
pct_bayes = mean(bayes) * 100.0

# Modal lag lengths by frequency
quarterly_lags = lags[freq .== 4]
monthly_lags = lags[freq .== 12]

function calc_mode(v)
    counts = Dict{eltype(v), Int}()
    for x in v
        counts[x] = get(counts, x, 0) + 1
    end
    return sort(collect(counts), by=x->x[2], rev=true)[1][1]
end

mode_q = Int(calc_mode(quarterly_lags))
mode_m = Int(calc_mode(monthly_lags))

println("===================================================================")
println("SECTION 5.1: LITERATURE SURVEY REPLICATION (81 PAPERS, 2015-2025)")
println("===================================================================")
@printf("Total Papers Surveyed: %d in top-6 economics journals\n", n_papers)
@printf("Modal lag length (quarterly): %d (Paper text: 4)\n", mode_q)
@printf("Modal lag length (monthly):   %d (Paper text: 12)\n", mode_m)
@printf("Lag Length / Frequency:\n")
@printf("   - Average across papers:   %.2f (Paper text: 0.96)\n", avg_lags_freq)
@printf("   - With Information Crit.:  %.2f (Paper text: 0.83)\n", ic_lags_freq)
@printf("Lag Length / Max Horizon:\n")
@printf("   - Average across papers:   %.2f%% (Paper text: 28%%)\n", avg_lags_maxhor * 100.0)
@printf("   - With Information Crit.:  %.2f%%\n", ic_lags_maxhor * 100.0)
@printf("Methodology Breakdown:\n")
@printf("   - Data-dependent lag selection (IC): %.1f%% (Paper text: 20%%)\n", pct_ic)
@printf("   - Bayesian shrinkage:                %.1f%% (Paper text: ~40%%)\n", pct_bayes)
println("===================================================================")

# Plot literature review histograms
Plots.default(fontfamily="sans-serif", dpi=300)

p_hist1 = histogram(lags_freq_rel, bins=10, normalize=:probability,
    color=:grey, alpha=0.7, label=false,
    title="Lag Length / Frequency",
    xlabel="Lag Length / Frequency", ylabel="Fraction")
vline!(p_hist1, [avg_lags_freq], color=:black, linewidth=2.5, label="Avg: $(round(avg_lags_freq, digits=2))")
vline!(p_hist1, [ic_lags_freq], color=:green, linestyle=:dash, linewidth=2.5, label="IC: $(round(ic_lags_freq, digits=2))")

p_hist2 = histogram(lags_maxhor_rel, bins=10, normalize=:probability,
    color=:grey, alpha=0.7, label=false,
    title="Lag Length / Max Horizon",
    xlabel="Lag Length / Max Horizon", ylabel="Fraction")
vline!(p_hist2, [avg_lags_maxhor], color=:black, linewidth=2.5, label="Avg: $(round(avg_lags_maxhor, digits=2))")
vline!(p_hist2, [ic_lags_maxhor], color=:green, linestyle=:dash, linewidth=2.5, label="IC: $(round(ic_lags_maxhor, digits=2))")

p_all = plot(p_hist1, p_hist2, layout=(1, 2), size=(880, 360), left_margin=8Plots.mm, bottom_margin=6Plots.mm)
fig_path = normpath(joinpath(@__DIR__, "..", "figures", "figure_lit_survey.png"))
savefig(p_all, fig_path)
println("Saved literature survey figure to: ", fig_path)
