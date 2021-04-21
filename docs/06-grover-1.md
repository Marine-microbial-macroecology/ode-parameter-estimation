# Turing-MCMC parameter estimation

```julia
using DataFrames
using CSV
using Plots
liefer = CSV.File("/Users/airwin/Dropbox/Julia/ode-parameter-estimation/liefer-growth-data.csv") |> DataFrame ; 
# liefer = CSV.File("liefer-growth-data.csv") |> DataFrame ; 
ss = filter( [:"Species", :"Replicate", :"Cell Density"] => (x,y,z) -> x == "Thalassiosira pseudonana" && y == "A" && !ismissing(z), liefer)
t = ss."Days in N-free Media"
R = ss.DIN
Q = ss.N
X = ss."Cell Density"
R = map(x -> ismissing(x) ? 0 : x, R)
RpgmL = R .* (10^3 * 14)
RpgmL[1] = RpgmL[2] + Q[2]*X[2] - Q[1]*X[1]
liefer[4,:"Dilution Factor"]
X[1] = X[1]*liefer[4,:"Dilution Factor"]
[t R Q X] # units: d, µmol/L, pg/cell, cells/mL
[t RpgmL Q X] # units: d, pg/mL, pg/cell, cells/mL
```

```
7×4 Matrix{Union{Missing, Real}}:
  0  -22687.9  2.6301    66851.2
  1   82028.8  3.19904       1.26367e5
  2       0    2.18655  218450.0
  3       0.0  1.57245       2.90062e5
  5       0    1.29513       3.52467e5
  7       0.0  1.11931       3.64933e5
 10       0    1.28601  316200.0
```





Look at the data

```julia
Plots.scatter(t, RpgmL .+ Q .* X, label = "Mass")
Plots.scatter!(t, RpgmL, label = "DIN")
Plots.scatter!(t, Q .* X, label = "Cell N")
```

![](data:image/png;base64,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)



Following example [tutorial](https://github.com/TuringLang/TuringTutorials/blob/master/10_diffeq.ipynb)
and [another source](https://turing.ml/dev/tutorials/10-bayesiandiffeq/)

```julia
using Turing, Distributions, DifferentialEquations 

# Import MCMCChain, Plots, and StatsPlots for visualizations and diagnostics.
using MCMCChains, Plots, StatsPlots

# Set a seed for reproducibility.
using Random
Random.seed!(14);
using Logging
Logging.disable_logging(Logging.Warn)
```

```
LogLevel(1001)
```





Define a model for the data

```julia
Turing.setadbackend(:forwarddiff)

@model function fitDroop(t, R, Q, X, logX)
    σ1 ~ InverseGamma(2, 3) # ~ is the tilde character
    σ2 ~ InverseGamma(2, 3) # ~ is the tilde character
    Km ~ truncated(Normal(100,10),0,200)
    Vmax ~ truncated(Normal(1.2,0.5),0,3)
    Qmin ~ truncated(Normal(1.0,0.5),0,3)
    muMax ~ truncated(Normal(1.0,0.5),0,3)

    p = [ Km, Vmax, Qmin, muMax]

    # must define the problem with numeric values first, then update with distributions
    prob1 = ODEProblem(droop!, [RpgmL1[1], Q1[1], X1[1]], (0.0, 10.0), [200.0, 1.0, 1.0, 1.0])
    prob = remake(prob1, p=p)  # modifies the original problem

    predicted = solve(prob, Rosenbrock23(), saveat=t)
    
    for j = 1:7
        Q[j] ~ Normal(predicted[j][2], σ1)
        logX[j] ~ Normal(log.(predicted[j][3]), σ2)
    end
end


@model function fitDroop1(t, R, Q, X, logX)
    σ1 ~ InverseGamma(2, 3) # ~ is the tilde character
    # σ2 ~ InverseGamma(2, 3) 
    R0 ~ Normal(300000, 1000)
    Q0 ~ truncated(Normal(3, 1), 0, 10)
    X0 ~ Normal(65000,1000)
    Km ~ truncated(Normal(100,10),0,200)
    Vmax ~ truncated(Normal(1.2,0.5),0,3)
    Qmin ~ truncated(Normal(1.0,0.5),0,3)
    muMax ~ truncated(Normal(1.0,0.5),0,3)

    p = [ Km, Vmax, Qmin, muMax]

    # must define the problem with numeric values first, then update with distributions
    prob1 = ODEProblem(droop!, [RpgmL1[1], Q1[1], X1[1]], (0.0, 10.0), [200.0, 1.0, 1.0, 1.0])
    prob = remake(prob1, u0=[R0, Q0, X0], p=p)  # modifies the original problem  # fails ****

    # prob = ODEProblem(droop!, [R0, Q0, X0], (0,10), p)
    # prob = ODEProblem(droop!, [R[1], Q[1], exp(X[1])], (0.0, 10.0), p)
    predicted = solve(prob, Rosenbrock23(), saveat=t)
    
    for j = 1:7
        Q[j] ~ Normal(predicted[j][2], σ1)
        # logX[i] ~ Normal(predicted[i][3], σ2)
    end
end
```

```
fitDroop1 (generic function with 1 method)
```





Create the model and use MCMC to fit it.


```julia
model = fitDroop(t, RpgmL1, Q1, X1, log.(X1))
chain2 = sample(model, NUTS(.65), MCMCThreads(), 100, 4, progress=false) # not enough iterations; demo only
```

<pre class="julia-error">
ERROR: UndefVarError: X1 not defined
</pre>




Median of posterior distribution

```julia
median(chain2[:muMax]), median(chain2[:Qmin]), median(chain2[:Km]), median(chain2[:Vmax])
```

<pre class="julia-error">
ERROR: UndefVarError: chain2 not defined
</pre>




Traceplots and distributions

```julia
Plots.plot(chain2)
```

<pre class="julia-error">
ERROR: UndefVarError: chain2 not defined
</pre>




Plots of data and solutions

```julia
chain_array = Array(chain2);

sol2 = solve(remake(prob, 
        p = [median(chain2[:Km]), median(chain2[:Qmin]), median(chain2[:Vmax]), median(chain2[:muMax])]), 
        Rosenbrock23()); 

pl = Plots.scatter(t, RpgmL1);
for k in 1:300
    resol = solve(remake(prob,p=chain_array[rand(1:size(chain_array)[1]), 1:4]),Rosenbrock23()) 
    # Note that due to a bug in AxisArray, the variables from the chain will be returned always in
    # the order it is stored in the array, not by the specified order in the call - :α, :β, :γ, :δ
    plot!(resol, vars=(0,1), alpha=0.3, color = "#BBBBBB", legend = false, ylims=(0, Inf))
end
plot!(sol2, vars=(0,1), alpha=1, color = "#BB0000", legend = false, ylims=(0, Inf))
display(pl)

pl = Plots.scatter(t, Q1);
for k in 1:300
    resol = solve(remake(prob,p=chain_array[rand(1:size(chain_array)[1]), 1:4]),Rosenbrock23()) 
    # Note that due to a bug in AxisArray, the variables from the chain will be returned always in
    # the order it is stored in the array, not by the specified order in the call - :α, :β, :γ, :δ
    plot!(resol, vars=(0,2), alpha=0.31, color = "#BBBBBB", legend = false)
end
plot!(sol2, vars=(0,2), alpha=1, color = "#BB0000", legend = false)
display(pl)

pl = Plots.scatter(t, log.(X1));
for k in 1:300
    resol = solve(remake(prob,p=chain_array[rand(1:size(chain_array)[1]), 1:4]),Rosenbrock23()) 
    # Note that due to a bug in AxisArray, the variables from the chain will be returned always in
    # the order it is stored in the array, not by the specified order in the call - :α, :β, :γ, :δ
    plot!(resol, vars=((t,x) -> (t, log.(x)), 0,3), alpha=0.3, color = "#BBBBBB", legend = false)
end
plot!(sol2, vars=((t,x) -> (t, log.(x)), 0,3), alpha=1, color = "#BB0000", legend = false)
display(pl)
```

<pre class="julia-error">
ERROR: UndefVarError: chain2 not defined
</pre>

