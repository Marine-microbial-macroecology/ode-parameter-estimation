# Parameter estimation from simulated data

## Pre-requisites

Julia packages used:

```@example Ex1
using DifferentialEquations, Plots
using Flux, DiffEqFlux, Optim, DiffEqSensitivity
import Statistics
using Turing, Distributions, DifferentialEquations 
using MCMCChains, StatsPlots
using Random
using Logging
```

## Simulating data

We use the same equations as before. This time we specify some times to sample the solution to obtain data to use in parameter estimation.


```@example Ex1
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

The initial conditions, parameters, and time-span for the solution must be specified.

```@example Ex1
u0 = [1.0, 1.0, 1.0]
p = [0.1, 2.0, 1.0, 0.8, 0.0, 0.0]
tspan = (0.0, 10.0)
tsteps = [0.1, 3.2, 4.5, 7.0, 9.1]
nothing # hide
```

Now we create and solve the ODE initial value problem. We will create two solutions,
one smoothly estimating the solution with interpolation on an interval and one with 
samples at a few discrete points.


```@example Ex1
prob = ODEProblem(droop!, u0, tspan, p)
sol1 = solve(prob, Tsit5())
sol2 = solve(prob, Rosenbrock23(), saveat = tsteps)
nothing # hide
```

We will plot the solution and the discrete samples.


```@example Ex1
Plots.plot(sol1)
Plots.scatter!(sol2)
```

We can convert the solution to a matrix and add some noise to the output.

```@example Ex1
data1 = Array(sol2) 
data1 = data1 .* (1 .+ 0.05*randn(size(data1)));
Plots.plot(sol1, label = ['R' 'Q' 'X'])
Plots.scatter!(sol2.t, data1', label = ['R' 'Q' 'X'])
```


## Optimization to find parameters

First we write a function to describe the difference between a trial solution and the
data points.  For the optimizer we use next, this should be a function of just the 
parameters to be adjusted.
Here my loss function is the difference between solution and data, scaled by standard deviation of each variable in data, squared and summed.

```@example Ex1
function loss(p)
    tspan = (0.0, 10.0)
    u0 = p[1:3]
    param = [ p[4:7] ; 0.0 ; 0.0 ] # force d = Rin = 0.0
    prob = ODEProblem(droop!, u0, tspan, param)
    sol2 = solve(prob, Rosenbrock23(), saveat = tsteps)
    data2 = Array(sol2)
    # loss = sum(((1.0 ./ Statistics.std(data1, dims=2)') * (data2 .- data1) ) .^ 2 )
    # loss = sum(((1.0 ./ (1 .+ Statistics.std(data1, dims=2)')) * (data2 .- data1) ) .^ 2 )
    loss = sum( (data2 .- data1) .^ 2 )
  return loss , sol2
end
```

(Note: I originally standardized the three variables as in the commented loss definition above, but 
this resulted in a worse solution, because R was so close to 0 and has a standard deviation of 10^{-7}. Perhaps a solution is to add a small amount to standard deviations to prevent this distortion. My first 
attempt, also commented out, did led to instabilities.)

The optimizer allows us to provide a callback function to show the loss score, or make a plot, at each iteration.  Here's an example callback function.

```@example Ex1
callback = function (p, l, pred)
  display(l)
  # plt = plot(pred, ylim = (0, 6))
  # display(plt)
  # Tell sciml_train to not halt the optimization. If return true, then
  # optimization stops.
  return false
end
```

We can test the loss function.

```@example Ex1
loss([1.0; 1.0; 1.0; p])
```

We are now ready to use the optimizer. ADAM is the gradient-search method.


```@example Ex1
result_ode = DiffEqFlux.sciml_train(loss, [1.0, 1.0, 1.0, 0.1, 2.0, 1.0, 1.0, 0, 0],
                                    ADAM(0.001),  # this argument needs to be small or instability
                                    cb = callback,
                                    maxiters = 200)
```

Now we use the parameters from this optimization to solve the differential equation
and plot the solution.

```@example Ex1
prob = ODEProblem(droop!, result_ode.u[1:3], tspan, result_ode.u[4:end])
sol3 = solve(prob, Rosenbrock23())
nothing # hide
```

We will plot the solution and the discrete samples.


```@example Ex1
Plots.plot(sol1, label = ["Original R" 'Q' 'X'])
Plots.plot!(sol3, lw = 2, label = ["New R" 'Q' 'X'])
Plots.scatter!(sol2.t, data1', label = ["Data R" 'Q' 'X'])
```

Compare parameters.

```@example Ex1
xaxis = ["R0" "Q0" "X0" "Km" "Vmax" "Qmin" "muMax" "d" "R0"]
estimate = result_ode.u'
original = [u0; p]'
error = (estimate .- original ) ./ original
p1 = Plots.scatter(xaxis, estimate, legend = false, title = "Parameters")
p1 = Plots.scatter!(xaxis, original, legend = false)
p2 = Plots.scatter(xaxis[1:7], error[1:7], legend = false, title = "Relative error")
plot(p1, p2, layout = (2,1))
```

## Bayesian MCMC parameter estimation

Set a random seed to make the result reproducible. Select an option for Turing package.

```@example Ex1
Random.seed!(14)
Turing.setadbackend(:forwarddiff)
nothing # hide
```

Define a model for the parameters, including priors, solution of the ODE, and comparison between data and the solution.

```@example Ex1
@model function fitDroop1(t, R, Q, X)
    σ1 ~ InverseGamma(2, 3)  # positive support; parameters α, β; mean β/(α-1), here 3
    σ2 ~ InverseGamma(2, 3) 
    R0 ~ truncated(Normal(1, 1), 0, 5) # might need strictly positive values?
    Q0 ~ truncated(Normal(1, 1), 0, 5)
    X0 ~ truncated(Normal(1, 1), 0, 5)
    Km ~ truncated(Normal(4.0, 2), 0, 5)
    Vmax ~ truncated(Normal(1.2, 2), 0, 5)
    Qmin ~ truncated(Normal(1.0, 2), 0, 5)
    muMax ~ truncated(Normal(1.0, 2), 0, 5)

    p = [ Km, Vmax, Qmin, muMax, 0.0, 0.0]

    # must define the problem with numeric values first, then update with distributions
    prob1 = ODEProblem(droop!, [R[1], Q[1], X[1]], (0.0, 10.0), [0.1, 1.0, 1.0, 1.0, 0.0, 0.0])
    prob = remake(prob1, u0=[R0, Q0, X0], p=p)  
    predicted = solve(prob, Rosenbrock23(), saveat=t)
    
    for j = 1:length(t)
        Q[j] ~ Normal(predicted[j][2], σ1)
        X[j] ~ Normal(predicted[j][3], σ2)
    end
end
nothing # hide
```

Create the model and simulate the chains.

```@example Ex1
t = sol2.t
R = [v[1] for v in sol2.u]
Q = [v[2] for v in sol2.u]
X = [v[3] for v in sol2.u]

model = fitDroop1(t, R, Q, X)
chain2 = sample(model, NUTS(.65), MCMCThreads(), 200, 4, progress=false)  # multi-threaded
# chain2 = mapreduce(c -> sample(model, NUTS(.75), 200), chainscat, 1:4) # single thread
nothing # hide
```

Extract data so we can plot trajectories from a selection of parameter values from the posterior distribution. Parameters come out of the chains in alphabetical order, so I resequence them to be in the order: initial conditons, parameter values as used in ODE function.


```@example Ex1
chain_array0 = Array(chain2);
chain_array = chain_array0[: ,[4, 2, 6, 1, 5, 3, 7] ]  # R Q X Km Vmax Qmin muMax
nothing # hide
```


```@example Ex1
N = size(chain_array)[1]
N2 = Int(N/2)
pR = Plots.scatter(t, R);
for k in 1:300
    pars = [ chain_array[rand(N2:N), :];  0.0; 0.0 ] # append d, Rin
    resol = solve(remake(prob, u0 = pars[1:3], p = pars[4:end]), Rosenbrock23()) 
    pR = plot!(resol, vars=(0,1), alpha=0.3, color = "#BBBBBB", legend = false)
end
pR = plot!(sol1, vars=(0,1), alpha=1, color = "#BB0000", legend = false, ylims=(0, Inf),
        yguide = "R")
pR = Plots.scatter!(t, R);
pQ = Plots.scatter(t, Q);
for k in 1:300
    pars = [ chain_array[rand(N2:N), :];  0.0; 0.0 ] # append d, Rin
    resol = solve(remake(prob, u0 = pars[1:3], p = pars[4:end]), Rosenbrock23()) 
    pQ = plot!(resol, vars=(0,2), alpha=0.3, color = "#BBBBBB", legend = false)
end
pQ = plot!(sol1, vars=(0,2), alpha=1, color = "#BB0000", legend = false, ylims=(0, Inf),
        yguide = "Q")
pQ = Plots.scatter!(t, Q);
pX = Plots.scatter(t, log.(X));
for k in 1:300
    pars = [ chain_array[rand(N2:N), :];  0.0; 0.0 ] # append d, Rin
    resol = solve(remake(prob, u0 = pars[1:3], p = pars[4:end]), Rosenbrock23()) 
    pX = plot!(resol, vars=((t,x) -> (t, log.(x)), 0, 3), alpha=0.3, color = "#BBBBBB", legend = false)
end
pX = plot!(sol1, vars=((t,x) -> (t, log.(x)), 0,3), alpha=1, color = "#BB0000", legend = false,
        yguide = "log X")
pX = Plots.scatter!(t, log.(X));
plot(pR, pQ, pX, layout = (3,1))
```

It looks like a few points leaves a lot of uncertainty about the trajectory!

