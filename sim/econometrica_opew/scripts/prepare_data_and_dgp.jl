# scripts/prepare_data_and_dgp.jl
# Loads Känzig (2021) data, cleans it, saves it to CSV, and estimates the VAR(18) population DGP

using CSV, DataFrames, LinearAlgebra, Statistics

data_dir = normpath(joinpath(@__DIR__, "..", "data"))
kaenzig_dir = joinpath(data_dir, "kaenzig")

println("Loading Känzig (2021) data...")
df_data = CSV.read(joinpath(kaenzig_dir, "OilDataM.csv"), DataFrame)
df_proxy = CSV.read(joinpath(kaenzig_dir, "OilSurprisesMLog.csv"), DataFrame, header=false)

poil = df_data.POIL[169:end]
cpi = df_data.CPI[169:end]
oilprod = df_data.OILPROD[169:end]
oilstocks = df_data.OILSTOCKS[169:end]
worldip = df_data.WORLDIP[169:end]
ip = df_data.IP[169:end]

# Variables transformation exactly as in Känzig (2021) and OPEW (2026):
# 1. log real oil price = log(POIL)*100 - log(CPI/100)*100
# 2. log world oil production = log(OILPROD)*100
# 3. log world oil inventories = log(OILSTOCKS)*100
# 4. log world industrial production = log(WORLDIP)*100
# 5. log US industrial production = log(IP)*100
# 6. log US CPI = log(CPI)*100
data = hcat(
    log.(poil) .* 100.0 .- log.(cpi ./ 100.0) .* 100.0,
    log.(oilprod) .* 100.0,
    log.(oilstocks) .* 100.0,
    log.(worldip) .* 100.0,
    log.(ip) .* 100.0,
    log.(cpi) .* 100.0
)

# Proxy: column 15 of OilSurprisesMLog
proxy = Float64.(df_proxy[!, 15])

# Combined 7-variable system: [proxy, data]
Y = hcat(proxy, data)
data_oil = Y .- mean(Y, dims=1)

# Save cleaned CSV for persistence
df_clean = DataFrame(
    oil_proxy = data_oil[:, 1],
    rpoil     = data_oil[:, 2],
    woprod    = data_oil[:, 3],
    woinv     = data_oil[:, 4],
    wip       = data_oil[:, 5],
    usip      = data_oil[:, 6],
    uscpi     = data_oil[:, 7]
)
clean_csv_path = joinpath(data_dir, "kanzig_cleaned_data.csv")
CSV.write(clean_csv_path, df_clean)
println("Saved cleaned Känzig dataset to: ", clean_csv_path)
println("Cleaned dataset dimensions: ", size(df_clean), " (528 months)")

# Estimate VAR(18) population DGP
p_pop = 18
T_raw, n_y = size(data_oil)
T_eff = T_raw - p_pop

Y_lhs = data_oil[(p_pop+1):end, :]
X_rhs = zeros(T_eff, n_y * p_pop)
for l in 1:p_pop
    X_rhs[:, ((l-1)*n_y + 1):(l*n_y)] = data_oil[(p_pop - l + 1):(end - l), :]
end

X_expand = hcat(X_rhs, ones(T_eff, 1))
betahat = Matrix((X_expand \ Y_lhs)') # n_y × (n_y * p_pop + 1)
B_est = betahat[:, 1:end-1]
const_est = betahat[:, end]
Resid = Y_lhs - X_expand * betahat'
Sigma_u = (Resid' * Resid) ./ (T_eff - size(betahat, 2))

# Companion matrix
A_c = zeros(n_y * p_pop, n_y * p_pop)
A_c[1:n_y, :] = B_est
if p_pop > 1
    A_c[(n_y+1):end, 1:(n_y*(p_pop-1))] = Matrix(1.0I, n_y*(p_pop-1), n_y*(p_pop-1))
end

evs = eigvals(A_c)
max_ev = maximum(abs.(evs))
println("VAR(18) estimated successfully.")
println("Max absolute eigenvalue of VAR(18): ", round(max_ev, digits=4))

