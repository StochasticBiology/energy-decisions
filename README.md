# energy-decisions
Simulations exploring the influence of energy supply on cellular decision-making

The simulation code here explores the behaviour of deterministic (ODE) and stochastic (Gillespie) models of a simple regulatory motif which exhibits different attractor basin structure under different parameterisations.

`grn-sim.c` run a collection of experiments, focussing on simulating the ODE system and outputting summary dynamics. It takes command line arguments specifying the experiment to run (to help parallelisation):

  * 0 -- Euler vs RK4
  * 1 -- param sweep with default model (takes another parameter 1-4 specifying the cooperativity structure of the model)
  * 2 -- param sweep with on-DNA dimerisation
  * 3 -- ATP x gamma3 influence
  * 4 -- ATP x cd2 influence
  * 5 -- different ICs
  * 6 -- different ATP influences
  * 7 -- bifurcation plots
  * 8 -- pairwise parameter influences (takes another parameter 0-7 specifying the pair of scalings to apply to the two parameters)

Split across cores, this ODE code should take at most a few hours to run. Experiment 1 will take a bit longer because one parameter set needs a rather smaller timestep for stability.

All the above produce CSV files describing the time behaviour of system state (protein levels, RNA levels, and so on). `grn-sim-plot.R` visualises these.

`stoch.c` simulates a stochastic version of the system using the Gillespie algorithm.

`stoch-scan.c` runs this stochastic version for the parameter sets in `param-scan.c`.

`hill.c` runs a simple simulation exploring the effects of different required oligomerisation states of a regulatory factor. `plot-hill.R` visualises the results for qualitative comparison with Hill function modelling.
 
