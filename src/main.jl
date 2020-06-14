using Pkg, Juno;
using CSV, DataFrames, Plots;
using Statistics, SparseArrays, LinearAlgebra;
cd(@__DIR__) #src

include("func_slp.jl")

mutable struct inputcom
  y        :: Array
  x        :: Array
  w        :: Array
  H_min    :: Int
  H_max    :: Int
  type     :: String
  r        :: Int8
  λ        :: Int8
end

df = CSV.read("../test/data.csv");
df = convert(Array{Float16}, df[:, 2:4]);
T, k = size(df);

P = 4; # Number of lags used in LP - controlled variable

# start LP at H_min=0 or 1 (H_min=1 if impose no contemporanous impact)
H_min = 1;
H_max = 20;

# Indicators for variables
ind_all = 1:k;
ind_response = 1; # Endogenous variable - response
ind_shock    = 3; # Endogenous variable related to the shock

ind_lp = "reg"; # Type of linear regression

ind_contempo = filter(x -> x ≠ ind_shock, ind_all);
#,(P+1):end
y  = df[:, ind_response]; # endogenous variable
x  = df[:, ind_shock]; # endoegnous variable related to the shock

# control variables (contemporaneous vars, lagged vars)
w  = [ df[:, ind_contempo]  matlag( df , P ) ];
w[.!isfinite.(w)] .= 0;
# Packaing everything as input
r = 2; #(r-1)=order of the limit polynomial (so r=2 implies the IR is shrunk towards a line )
λ = 100; # value of penalty

testing₀ = inputcom(y, x, w, H_min, H_max, "reg", 3,4);
testing₁ = inputcom(y, x, w, H_min, H_max, "smooth", r, λ);

lp₀ = smoothlocalprojection(testing₀);
#Juno.@enter smoothlocalprojection(testing₀);
lp₁ = smoothlocalprojection(testing₁);

plot(lp₀.IR)
plot!(lp₁.IR)

# Print responses in simple text format
notes = ["C4", "D4", "E4", "F4"]
outfile = "notes.txt"
open(outfile, "w") do f
  for i in notes
    println(f, i)
  end
end # the file f is automatically closed after this block finishes
