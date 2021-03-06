# Posterior distribution of parameter 

Estimate parameters and initial conditions for the Droop-Grover system using data from Liefer et al (2019)
woth Turing and MCMC to get posterior distributions on the parameters.

Start by preparing the data the same way we did for the point-estimate example

```julia
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

Tp = filter([:Species, :cells ] => (x,y) -> x == "Thalassiosira pseudonana" && !ismissing(y), liefer)[:, [:days, :Replicate, :DIN_pgml, :N, :cells, :mass]]
Tw = filter([:Species, :cells ] => (x,y) -> x == "Thalassiosira weissflogii" && !ismissing(y), liefer)[:, [:days, :Replicate, :DIN_pgml, :N, :cells, :mass]]
Ot = filter([:Species, :cells ] => (x,y) -> x == "Ostreococcus tauri" && !ismissing(y), liefer)[:, [:days, :Replicate, :DIN_pgml, :N, :cells, :mass]]
Ms = filter([:Species, :cells ] => (x,y) -> x == "Micromonas sp." && !ismissing(y), liefer)[:, [:days, :Replicate, :DIN_pgml, :N, :cells, :mass]]

# define ODE
function droop!(du, u, p, t)
  R, Q, X = u
  Km, Vmax, Qmin, muMax = p # , d, R0 = p
  d = 0.0
  R0 = 0.0
  rho = Vmax * R / (Km + R)
  mu = muMax * (1 - Qmin/Q)
  du[1] = dRdt = d*(R0 - R) - rho*X
  du[2] = dQdt = rho - mu*Q
  du[3] = dXdt = (mu - d)*X
end
```




Create model to fit to data.

First simplified using one replicate.

```julia
@model function fitDroop0(DF)
    data0a = identity.(Array(filter(:Replicate => x -> x == "A", DF)[:, [:days, :DIN_pgml, :N, :cells]])') 
    t_a = data0a[1,:]
    R_a = data0a[2,:]
    Q_a = data0a[3,:]
    X_a = data0a[4,:]
    tmax = maximum(t_a)
    
    ??1 ~ InverseGamma(2, 3)
    ??2 ~ InverseGamma(2, 3) 
    ??3 ~ InverseGamma(2, 3) 
    R0a ~ truncated(Normal(R_a[1], std(R_a)), 0.0, 10.0*maximum(R_a))
    Q0a ~ truncated(Normal(Q_a[1], std(Q_a)), 0.0, 10.0*maximum(Q_a))
    X0a ~ Normal(X_a[1], std(X_a))

    Km ~ truncated(Normal(R_a[1]/5.0,10), 0.0, maximum(R_a))
    Vmax ~ truncated(Normal(1.2, 1.0), 0.0, 3)
    Qmin ~ truncated(Normal(Q_a[1]/10.0, 0.5), 0.0, maximum(Q_a))
    muMax ~ truncated(Normal(1.0,0.5), 0.0, 3)

    p = [ Km, Vmax, Qmin, muMax]

    # must define the problem with numeric values first, then update with distributions
    prob1 = ODEProblem(droop!, [R_a[1], Q_a[1], X_a[1]], (0.0, tmax), [200.0, 1.0, 1.0, 1.0])
    probA = remake(prob1, u0=[R0a, Q0a, X0a], p=p) 

    predictedA = solve(probA, Rosenbrock23(), saveat=t_a)
    
    for j = 1:length(t_a)
        # R_a[j] ~ Normal(predictedA[j][1], ??1)
        Q_a[j] ~ Normal(predictedA[j][2], ??2)
        # X_a[j] ~ Normal(predictedA[j][3], ??3)
    end
end
```




Now using all three replicates.

```julia
@model function fitDroop1(DF)
    data0a = identity.(Array(filter(:Replicate => x -> x == "A", DF)[:, [:days, :DIN_pgml, :N, :cells]])') 
    data0b = identity.(Array(filter(:Replicate => x -> x == "B", DF)[:, [:days, :DIN_pgml, :N, :cells]])')
    data0c = identity.(Array(filter(:Replicate => x -> x == "C", DF)[:, [:days, :DIN_pgml, :N, :cells]])')
    t_a = data0a[1,:]
    R_a = data0a[2,:]
    Q_a = data0a[3,:]
    X_a = data0a[4,:]
    #logX_a = log.(X_a)
    t_b = data0b[1,:]
    R_b = data0b[2,:]
    Q_b = data0b[3,:]
    X_b = data0b[4,:]
    #logX_b = log.(X_b)
    t_c = data0c[1,:]
    R_c = data0c[2,:]
    Q_c = data0c[3,:]
    X_c = data0c[4,:]
    #logX_c = log.(X_c)
    tmax = maximum([t_a t_b t_c])
    
    ??1 ~ InverseGamma(2, 3)
    ??2 ~ InverseGamma(2, 3) 
    ??3 ~ InverseGamma(2, 3) 
    R0a ~ truncated(Normal(R_a[1], std(R_a)), 0.0, 10.0*maximum(R_a))
    Q0a ~ truncated(Normal(Q_a[1], std(Q_a)), 0.0, 10.0*maximum(Q_a))
    X0a ~ Normal(X_a[1], std(X_a))
    R0b ~ truncated(Normal(R_b[1], std(R_b)), 0.0, 10.0*maximum(R_b))
    Q0b ~ truncated(Normal(Q_b[1], std(Q_b)), 0.0, 10.0*maximum(Q_b))
    X0b ~ Normal(X_b[1], std(X_b))
    R0c ~ truncated(Normal(R_c[1], std(R_c)), 0.0, 10.0*maximum(R_c))
    Q0c ~ truncated(Normal(Q_c[1], std(Q_c)), 0.0, 10.0*maximum(Q_c))
    X0c ~ Normal(X_c[1], std(X_c))

    Km ~ truncated(Normal(R_a[1]/5.0,10), 0.0, maximum(R_a))
    Vmax ~ truncated(Normal(1.2, 1.0), 0.0, 3)
    Qmin ~ truncated(Normal(Q_a[1]/10.0, 0.5), 0.0, maximum(Q_a))
    muMax ~ truncated(Normal(1.0,0.5), 0.0, 3)

    p = [ Km, Vmax, Qmin, muMax]

    # must define the problem with numeric values first, then update with distributions
    prob1 = ODEProblem(droop!, [R_a[1], Q_a[1], X_a[1]], (0.0, tmax), [200.0, 1.0, 1.0, 1.0])
    probA = remake(prob1, u0=[R0a, Q0a, X0a], p=p) 
    probB = remake(prob1, u0=[R0b, Q0b, X0b], p=p) 
    probC = remake(prob1, u0=[R0c, Q0c, X0c], p=p) 

    predictedA = solve(probA, Rosenbrock23(), saveat=t_a)
    predictedB = solve(probB, Rosenbrock23(), saveat=t_b)
    predictedC = solve(probC, Rosenbrock23(), saveat=t_c)
    
    for j = 1:length(t_a)
        R_a[j] ~ Normal(predictedA[j][1], ??1)
        Q_a[j] ~ Normal(predictedA[j][2], ??2)
        X_a[j] ~ Normal(predictedA[j][3], ??3)
    end
    for j = 1:length(t_b)
        R_b[j] ~ Normal(predictedB[j][1], ??1)
        Q_b[j] ~ Normal(predictedB[j][2], ??2)
        X_b[j] ~ Normal(predictedB[j][3], ??3)
    end
    for j = 1:length(t_c)
        R_c[j] ~ Normal(predictedC[j][1], ??1)
        Q_c[j] ~ Normal(predictedC[j][2], ??2)
        X_c[j] ~ Normal(predictedC[j][3], ??3)
    end
end
```




## Thalasiosira weisflogii

```julia
model0 = fitDroop0(Tw)
model1 = fitDroop1(Tw)
```




Fit model

```julia
chain0 = sample(model0, NUTS(0.65), 10)  
chain1 = sample(model1, NUTS(0.65), 100)  
# chain = mapreduce(c -> sample(model, NUTS(.65), 1000), chainscat, 1:4) # not multithreaded
# chain = sample(model, NUTS(.65), MCMCThreads(), 100, 4, progress=false) # multithreaded
```




Fitting the chains gives an error: ERROR: TypeError: in typeassert, expected Float64, got a value of type ForwardDiff.Dual{Nothing, Float64, 9}

to be resolved later on....


show results

```julia
chain0
```

<pre class="julia-error">
ERROR: UndefVarError: chain0 not defined
</pre>


```julia
chain1
```

<pre class="julia-error">
ERROR: UndefVarError: chain1 not defined
</pre>


```julia
plot(chain0)
```

<pre class="julia-error">
ERROR: UndefVarError: chain0 not defined
</pre>

