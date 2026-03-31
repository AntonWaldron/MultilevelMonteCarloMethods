% This is the test function using the level l estimator for pi.
% This function produces the convergence plots and results for the pi estimator used in Chapter 2 of the FYP. 

function mlmcpi_test

close all; clear all;

addpath('..');

N0 = 10000;
Lmin = 4;
Lmax = 8;

fprintf(1, '\n ---- Pi Estimator ---\n');
N = 1000;
L = 6;
Eps    = [ 0.005 0.01 0.02 0.05 0.1 ];

filename = 'MLMC_pi_Term2_M=8_n=100';

fp = fopen([filename '.txt'],'w');
mlmc_test(@pi_l, N,L, N0, Eps, Lmin, Lmax, fp)
fclose(fp);

nvert = 3;
mlmc_plot(filename, nvert);

if(nvert==1)
    figure(1)
    print('-dpdf',[filename 'a.pdf'])
    figure(2)
    print('-dpdf',[filename 'b.pdf'])
else
    print('-dpdf',[filename '.pdf'])
end

filename = [filename '_n100'];
fp = fopen([filename '.txt'], 'w');
% Now do 100 MLMC calcs in parallel
mlmc_test_100(@pi_l, pi, N0, Eps, Lmin, Lmax , fp);
fclose(fp);
mlmc_plot_100(filename);
print('-dpdf', [filename '.pdf'])
end
