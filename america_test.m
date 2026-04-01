% This function, combined with the america_multilevel.m and mlmc functions produce the convergence and result plots used in the American options section of the FYP. 

function america_test


% Out of memory for L > 4

close all; clear all;

addpath('..');

N0    = 5000;   % initial samples on coarse levels
Lmin  = 2;      % minimum refinement level
Lmax  = 5;      % maximum refinement level

% Adjust these maybe
fprintf(1,'\n ---- American put ---- \n');
N      = 2500;       % samples for convergence tests
L      = 5;             % levels for convergence tests 
Eps    = [ 0.02 0.05 0.1 0.2 ];

filename = 'America_Multilevel_M=4';
fp = fopen([filename '.txt'],'w');
mlmc_test_2(@america_multilevel, N,L, N0,Eps,Lmin,Lmax, fp);
fclose(fp);

% % Changed to pdf from eps
nvert = 4;
mlmc_plot_2(filename, nvert);

if(nvert==1)
    figure(1)
    print('-dpdf',[filename 'a.pdf'])
    figure(2)
    print('-dpdf',[filename 'b.pdf'])
else
   % print('-dpdf',[filename '.pdf'])
end

%
% now do 100 MLMC calcs in parallel
%
  filename = [filename 'fixed_n'];
  fp = fopen([filename '.txt'],'w');
  mlmc_test_100(@america_multilevel, NaN, N0,Eps,Lmin,Lmax, fp);
  fclose(fp);

%
% plot results
% 
  mlmc_plot_100(filename);
  print('-dpdf',[filename '.pdf'])
