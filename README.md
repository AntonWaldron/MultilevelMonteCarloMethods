# MultilevelMonteCarloMethods

This is a code repository containing all of the MATLAB M-functions used to produce results for my Final Year Project, titled 'Multilevel Monte Carlo Methods for Options Pricing'. 

The M-functions opre.m, mcqmc06.m, mlmc.m, mlmc_plot.m, mlmc_test.m, mlmc_test_100.m, are all downloaded from Giles' webpage https://people.maths.ox.ac.uk/gilesm/mlmc/. 

The M-functions america_multilevel.m and america_test.m were written in order to price American options using the Long-Staff and Schwarz method in a multilevel context.

The M-function calloptionspricing.m  reads in a datset of call options and produces estimates for the MLMC and Standard cost needed for a RMSE target of 0.01 for European options.

The BPATH files produce the Brownian Motion plots used in the report. 
