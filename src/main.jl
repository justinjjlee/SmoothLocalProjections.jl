using Pkg, Juno;
using CSV, DataFrames, Plots;
using Statistics, Distributions, SparseArrays, LinearAlgebra;
cd(@__DIR__) #src

include("func_slp.jl");


df = CSV.read("../test/data.csv");
df = convert(Array{Float64}, df[:, 2:4]);

## Parameterization
P = 4; # Number of lags used in LP - controlled variable

# start LP at H_min=0 or 1 (H_min=1 if impose no contemporanous impact)
H_max = 20;
ind_response = 1; # Endogenous variable - response
ind_shock    = 3; # Endogenous variable related to the shock

# Packaing everything as input
r = 2; #(r-1)=order of the limit polynomial
# NOTE: (so r=2 implies the IR is shrunk towards a line )
λ = 100; 

indx     = [ind_response ind_shock H_max P];

inicjał₀ = initalz(df, indx, "reg", [r 0]);
inicjał₁ = initalz(df, indx, "smooth", [r λ]);
# Obliczanie projekcja lokalna, in order of,
#  (1) Local projection a la Jordá (2005)
#  (2) Local projection with basic parameteri tested
#  (3) Local projection with optimal parameter validated
lp₀, lp₁ = slp(inicjał₀), slp(inicjał₁);

# groß - evaluation using validation exam,
# Cross-valudation for value of optimal value of λ → λₒ
λₒ, resvec = slpᵥ(lp₁);
plt = plot(resvec[:,1], resvec[:,2],
    title = "walidacja krzyżowa dla optymalny λ",
    xlabel = "λ",
    ylabel = "MSE",
    label = "oszacowanie")
plot!((λₒ .* ones(2, 1)),
      [minimum(resvec[:,2]) minimum(resvec[:,2])+((maximum(resvec[:,2])-minimum(resvec[:,2]))/2)]',
      label = "optymalny"
      )
# optymalny
inicjał₂ = initalz(df, indx, "smooth", [r λₒ]);
lp₂ = slp(inicjał₂);
lp₂_ci = slp_ci(lp₂);

# rezultat ============================================================================================
plot(lp₀.IR, xlabel = "Time since stimulus/impact", ylabel = "Response", label = "Jordá (2005)")
plot!(lp₁.IR, label = "SLP: λ = $(λ)", color = "purple")
plot!(zeros(length(lp₁.IR)), label = false, color = "black", line = :dot)
plot!(lp₂.IR, label = "SLP: λ optymalny = $(λₒ)", color = "red")
plot!(lp₂_ci, line = :dash, color = "red", legend = false)

# fantazyjny ==========================================================================================
function plt_anime(ir₁, ir₂, iter)
  plot!(ir₁.IR[1:iter], label = "SLP: λ = $(λ)", color = "purple")
  plot!(ir₂.IR[1:iter], label = "SLP: λ optymalny = $(λₒ)", color = "red")
end

plot(lp₀.IR, xlabel = "Time since stimulus/impact", ylabel = "Response", label = "Jordá (2005)")
plot!(zeros(length(lp₁.IR)), label = false, color = "black", line = :dot)
plot!(lp₂_ci, line = :dash, color = "red", legend = false)
fantazyjny = @animate for iter ∈ 1:length(lp₁.IR)
  plt_anime(lp₁, lp₂, iter)
end

gif(fantazyjny, "example.gif", fps = 5)