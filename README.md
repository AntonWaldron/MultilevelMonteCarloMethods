# MultilevelMonteCarloMethods

This is a code repository containing all of the MATLAB M-functions used to produce results for my Final Year Project titled 'Multilevel Monte Carlo Methods for Options Pricing'. 

The M-functions opre.m, mcqmc06.m, mlmc.m, mlmc_plot.m, mlmc_test.m, mlmc_test_100.m, are all downloaded from Giles' webpage https://people.maths.ox.ac.uk/gilesm/mlmc/. They are used, in conjuction with edits to the source code from Giles, to produce the majority of the plots and results the report. 

The M-functions pi_l and mlmcpi_test provide a nice introduction to MLMC in application; the first shows how to adapt a standard Monte Carlo approach to produce a level l estimator and the second calls the MLMC driver functions to produce the convergence and result plots used in the final report. 

The M-functions america_multilevel.m and america_test.m were written to price American options using the Long-Staff and Schwarz method in a multilevel context.

The M-function calloptionspricing.m  reads in a dataset of call options and produces estimates for the MLMC and Standard cost needed for an RMSE target of 0.01 for European options.

The BPATH files are adapted from D. Higham's 'An Algorithmic Introduction to Numerical Simulation of Stochastic Differential Equations' paper and are used to produce the Brownian motion, Brownian Bridging, and the illustration of the multilevel coupling for Brownian Motions plots used in the report. 
