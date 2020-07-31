using Pkg, Juno;
using CSV, DataFrames, Plots;
using Statistics, SparseArrays, LinearAlgebra;
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
λ = 100; # value of penalty

indx     = [ind_response ind_shock H_max P];
param    = [r λ];

testing₀ = initalz(df, indx, "reg", param);
testing₁ = initalz(df, indx, "smooth", param);
## Point estimation of the local projection
lp₀ = slp(testing₀);
#Juno.@enter smoothlocalprojection(testing₀);
lp₁ = slp(testing₁);
# COmpare plots
plot(lp₀.IR)
plot!(lp₁.IR)

## Cross-validation of optimal λ
# Cross-valudation for value of optimal value of λ
λₒ, resvec = slpᵥ(lp₁);
#Juno.@enter slpᵥ(lp₁)
plot(resvec[:,1], resvec[:,2])



# Print responses in simple text format
notes = ["C4", "D4", "E4", "F4"]
outfile = "notes.txt"
open(outfile, "w") do f
  for i in notes
    println(f, i)
  end
end # the file f is automatically closed after this block finishes
