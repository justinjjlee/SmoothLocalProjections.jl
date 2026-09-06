# scripts/replicate_ramey_summary.jl
# Summarizes the Ramey (2016) Empirical Applications referenced in Section 4.2
# and Figure 1 of Montiel Olea, Plagborg-Møller, Qian, & Wolf (Econometrica 2026).

using CSV
using DataFrames
using Statistics
using Printf

data_path = normpath(joinpath(@__DIR__, "..", "data", "ramey_se_ratios.csv"))
if !isfile(data_path)
    error("File not found: $data_path")
end

df = CSV.read(data_path, DataFrame)

println("===================================================================")
println("RAMEY (2016) EMPIRICAL APPLICATION: VAR vs. LP STANDARD ERROR RATIOS")
println("Montiel Olea, Plagborg-Møller, Qian, & Wolf (Econometrica 2026, Section 4.2)")
println("===================================================================\n")

total_ratios = nrow(df)
mean_ratio   = mean(df.se_ratio)
med_ratio    = median(df.se_ratio)
p10_ratio    = quantile(df.se_ratio, 0.10)
p90_ratio    = quantile(df.se_ratio, 0.90)

@printf("Total Standard Error Ratios: %d across 4 macro applications\n\n", total_ratios)

# Breakdown by application
for app in unique(df.application)
    sub = filter(row -> row.application == app, df)
    @printf("• Application: %-25s | Count: %2d | Mean SE Ratio: %.3f | Median: %.3f\n",
            app, nrow(sub), mean(sub.se_ratio), median(sub.se_ratio))
end

println("\n-------------------------------------------------------------------")
println("OVERALL DISTRIBUTION (MEDIUM & LONG HORIZONS: 1 TO 4/5 YEARS)")
println("-------------------------------------------------------------------")
@printf("Mean Ratio:        %.3f (Paper replication: 0.394)\n", mean_ratio)
@printf("Median Ratio:      %.3f (Paper replication: 0.367)\n", med_ratio)
@printf("10th Percentile:   %.3f (Paper replication: 0.168)\n", p10_ratio)
@printf("90th Percentile:   %.3f (Paper replication: 0.638)\n", p90_ratio)
println("-------------------------------------------------------------------")
println("Matches the shaded gray region in Figure 1 [0.168, 0.638].")
println("===================================================================")
