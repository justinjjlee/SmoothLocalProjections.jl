using Pkg;
using CSV, DataFrames, Plots;
using Statistics, Distributions, SparseArrays, LinearAlgebra;
using Test

cd(@__DIR__) #src

include("../src/functions.jl")

@testset "SmoothLocalProjections.jl" begin
    df = CSV.read("data.csv", DataFrame);
    df = convert(Array{Float16}, df[:, 2:4]);
    # .....................................................
    # Parameters to be used
    T, k = size(df);
    P = 4; # Number of lags used in LP - controlled variable
    # start LP at H_min=0 or 1 (H_min=1 if impose no contemporanous impact)
    H_max = 20;
    ind_response = 1; # Endogenous variable - response
    ind_shock    = 3; # Endogenous variable related to the shock
    # Packaing everything as input
    r = 2; #(r-1)=order of the limit polynomial
    # NOTE: (so r=2 implies the IR is shrunk towards a line )
    λ = 100; 
    # Combine parameters
    indx     = [ind_response ind_shock H_max P];
    # .....................................................
    # Test the data loads
    @test typeof(df) <: Matrix
    @test any(indx .== NaN) == false

    # Run through the sample exercise
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

    # Optimal parameter
    inicjał₂ = initalz(df, indx, "smooth", [r λₒ]);
    lp₂ = slp(inicjał₂);
    lp₂_ci = slp_ci(lp₂);


    @test length(lp₀.IR) == (H_max - 0 + 1)
    @test any(lp₀.IR .== NaN) == false
    @test any(lp₁.IR .== NaN) == false
end