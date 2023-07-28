# energy-decisions
Simulations exploring the influence of energy supply on cellular decision-making

The simulation code here explores the behaviour of deterministic (ODE) and stochastic (Gillespie) models of a simple regulatory motif which exhibits different attractor basin structure under different parameterisations.

`param-scan.c` sweeps parameter space of the system, plotting steady-state behaviour as a function of initial conditions for different rate parameters and degrees of cooperativity.

`zoom-scan.c` explores a more specific set of parameters corresponding to rates determining attractor basin structure and energy availability.

`time-series.c` outputs explicit solutions as time series for a representative set of parameters.

`stoch.c` simulates a stochastic version of the system using the Gillespie algorithm.

All the above produce CSV files describing the time behaviour of system state (protein levels, RNA levels, and so on). `plot-vis.R` takes output from all the above and produces visualisations of the system behaviour.

`hill.c` runs a simple simulation exploring the effects of different required oligomerisation states of a regulatory factor. `plot-hill.R` visualises the results for qualitative comparison with Hill function modelling.
