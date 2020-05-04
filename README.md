# julia-SmoothLocalProjections.jl
Implementation of Smooth Local Projections (SLP) based on [Barnichon and Brownlees (2019)](https://www.mitpressjournals.org/doi/abs/10.1162/rest_a_00778) - "Impulse Response Estimation by Smooth Local Projections." Original method of Local Projections is first introduced by [Òscar Jordà (2005)](https://www.aeaweb.org/articles?id=10.1257/0002828053828518)

The code was translated from MATLAB code published from the replication file made available by [C. Brownlees](https://github.com/ctbrownlees/MATLAB-package-lproj).

```julia
import Pkg;
Pkg.update();

using CSV, Statistics, LinearAlgebra: I;
cd(@__DIR__) #src location

include("func_slp.jl")
```

Input argument is constructed using a mutable object of following,

```julia
mutable struct InputPackage
  y        :: Array        # Defined endogenous variable for response.;
  x        :: Array        # Defined endogenous variable for generating shocks.;
  w        :: Array        # Defined endogenous variable - contemporaneous and lagged structure.;
  H_min    :: Int          # Horizon sought - start.;
  H_max    :: Int          # Horizon sought - end.;
  type     :: String       # Type of regression: 'reg' to be regular, smooth otherwise.;
  r        <: Float64      # Order of limit polynomial: 1 for smooth, 2 for linear pattern.;
  λ        <: Float64      # Value of penalty: Barnichon and Brwonlees (2019) uses 100.;
end
```
