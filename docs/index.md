# Parameter estimation for ordinary differential equations

In these notes I will demonstrate some techniques for estimating parameters and initial 
conditions for ordinary differential equations. The application will be the Droop-Grover 
model of phytoplankton growth on one limiting nutrient with growth rate determined by internal nutrient reserves.

Examples are broken into the following sections

* [Section 1](02-eqns.html): The differential equations and their numerical solution
* [Section 2](03-simulations.html): Parameter estimation from simulated data using both "Flux" optimization and "Turing" model with MCMC
* [Section 3](04-read-data.html): Gather some observational data from four lab cultures undergoing nitrogen starvation
* [Section 4](05-data-point-estimate.html): SciML optimization to estimate parameters and initial conditions for these data
* [Section 5](06-grover-1.html): Turing MCMCChains estimation of posterior distributions for parameters (early version)
* [Section 5](06-posterior-distribution.html): (not working yet) Turing-MCMCChains methods to estimate distributions of parameters and predicted trajectories

A related effort is to discover the differential equations that drive a system by sparse regression. That is, provide many possible terms for the right hand side of $x' = f(t,x)$ and 
use LASSO-type shrinkage to eliminate most of the terms.

* [Section 2.1](07-data-drive-de.html): An example using the Lorenz equations
* [Section 2.2](xxx): not written yet. Droop-Grover phytoplankton growth Model
* [Section 2.3](xxx): An example using culture data (not written yet)

Here are a couple of other ideas related to growth of phytoplankton cultures.

* [Section 3.1](20-growth-gam.html): Estimating growth rate in a non-steady-state culture using a model-free (GAM) approach
* [Section 3.2](21-growth-glmnet.html): Attempt to discover the macromolecular ratios that best predict phytoplankton growth under N starvation using sparse regression (GLMnet/LASSO)

Some other Julia explorations

* [Section X.1](08-probabilistic-solutions.html): Explore trajectories of numerical soltutions of ODEs. (A quick example to explore more fully later.)




## Questions

I have lots of Julia questions!

* What is a good way to refer to files relative to the path of the current file or project? I've used absolute paths in these files, which seems like a terrible idea
* What optimizer am I using with scikit_optim, ADAM ? Should I be doing something else?
* What is a good way to cope with defining an error/RMSE function with: one variable very close to 0 for most but not all observations, one variable lognormal, and one variable normally distributed?

And some other questions!

* Math shows up fine in HTML, but is not formatted in md / github pages
