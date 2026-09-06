# scripts/test_vma.jl
using CSV, DataFrames, LinearAlgebra, Statistics, Distributions

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

println("VMA(1) matches C_pop: ", norm(VMA_IRF[1, :, :] - C_pop))
println("VMA(2) norm: ", norm(VMA_IRF[2, :, :]))
println("VMA(13) norm (omitted lag 13 effect): ", norm(VMA_IRF[13, :, :]))
