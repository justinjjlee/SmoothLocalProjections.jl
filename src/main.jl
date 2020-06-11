using Pkg, Juno;
using CSV, Statistics, LinearAlgebra: I;
cd(@__DIR__) #src

include("func_slp.jl")

mutable struct InputPackage
  y        :: Array
  x        :: Array
  w        :: Array
  H_min    :: Int
  H_max    :: Int
  type     :: String
  r        <: Float64
  λ        <: Float64
end

df = CSV.read("data.csv");
df = convert(Array{Float64}, df[:, 2:4]);
T, k = size(df);

P = 4; # Number of lags used in LP - controlled variable

# start LP at H_min=0 or 1 (H_min=1 if impose no contemporanous impact)
H_min = 1;
H_max = 20;

# Indicators for variables
ind_all = 1:k
ind_response = 1; # Endogenous variable - response
ind_shock    = 3; # Endogenous variable related to the shock

ind_lp = "reg"; # Type of linear regression

ind_contempo = filter(x -> x ≠ ind_shock, ind_all);

y  = df[(P+1):end, ind_response]; # endogenous variable
x  = df[(P+1):end, ind_shock]; # endoegnous variable related to the shock

# control variables (contemporaneous vars, lagged vars)
w  = [ df[(P+1):end, ind_contempo]  matlag( df , P ) ];

# Packaing everything as input

testing = InputPackage(y, x, w, H_min, H_max,ind_lp, 3,4)

# Edit by 11/13/2019



w( ~isfinite(w) ) = 0;

lp = locproj(y,x,w,H_min,H_max,ind_lp); % IR from (standard) Local Projection

r = 2; %(r-1)=order of the limit polynomial (so r=2 implies the IR is shrunk towards a line )
lambda = 100; % value of penalty
