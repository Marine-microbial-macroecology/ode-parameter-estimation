# Posterior distribution of parameter 

Of parameters an initial conditions for the Droop-Grover system using data from Liefer et al (2019).

Start by preparing the data the same way we did for the point-estimate example

```@setup Ex1   # hide
# packages
using DataFrames, CSV, StatsPlots, Query
using DifferentialEquations, Plots
using Flux, DiffEqFlux, Optim, DiffEqSensitivity
import Statistics
using LinearAlgebra

# tools for probabilistic modelling
using Turing, Distributions, DifferentialEquations 

# Import MCMCChain, Plots, and StatsPlots for visualizations and diagnostics.
using MCMCChains, Plots, StatsPlots

# Set a seed for reproducibility.
using Random
Random.seed!(14);
using Logging
Logging.disable_logging(Logging.Warn)

# read data
liefer = CSV.File("liefer-growth-data.csv") |> DataFrame ; 
rename!(liefer, [2 => :days, 5 => :cell_density, 7 => :dilution_factor])
replace!(liefer.dilution_factor, missing => 1.0)
replace!(liefer.DIN, missing => 0.0)
transform!(liefer, :DIN => (x -> x .* (10^3 * 14)) => :DIN_pgml)
gdf = groupby(liefer, [:Species, :Replicate])
transform!(gdf, :dilution_factor => cumprod)
transform!(gdf, :dilution_factor_cumprod => ( (x) -> x ./ minimum(x)) => :dilution_factor_min)
transform!(gdf, [:cell_density, :dilution_factor_min] => ( ./ ) => :cells)
transform!(liefer, [:days, :DIN_pgml] => ByRow((x,y) -> x === 0 ? missing : y) => :DIN_pgml_corr)
transform!(liefer, [:DIN_pgml_corr, :N, :cells] => ( (R,Q,X) -> R .+ Q .* X ) => :mass)

Tp = filter([:Species, :cell_density ] => (x,y) -> x == "Thalassiosira pseudonana" && !ismissing(y), liefer)[:, [:days, :Replicate, :DIN_pgml, :N, :cell_density, :mass]]
Tw = filter([:Species, :cell_density ] => (x,y) -> x == "Thalassiosira weissflogii" && !ismissing(y), liefer)[:, [:days, :Replicate, :DIN_pgml, :N, :cell_density, :mass]]
Ot = filter([:Species, :cell_density ] => (x,y) -> x == "Ostreococcus tauri" && !ismissing(y), liefer)[:, [:days, :Replicate, :DIN_pgml, :N, :cell_density, :mass]]
Ms = filter([:Species, :cell_density ] => (x,y) -> x == "Micromonas sp." && !ismissing(y), liefer)[:, [:days, :Replicate, :DIN_pgml, :N, :cell_density, :mass]]

# define ODE
function droop!(du, u, p, t)
  R, Q, X = u
  Km, Vmax, Qmin, muMax, d, R0 = p
  rho = Vmax * R / (Km + R)
  mu = muMax * (1 - Qmin/Q)
  du[1] = dRdt = d*(R0 - R) - rho*X
  du[2] = dQdt = rho - mu*Q
  du[3] = dXdt = (mu - d)*X
end

nothing # hide
```
