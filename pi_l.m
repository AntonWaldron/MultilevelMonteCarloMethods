%% This is a level l estimator function for pi_l used in the introductory section of the FYP. 

function [sums, cost] = pi_l(l, N)

%mlmc_l = function for level l estimator 

% [sums, cost] = mlmc_fn(l,N, varargin)     low-level routine
%
% inputs:  l = level
%          N = number of samples
%          varargin = optional additional user variables
%
% output: sums(1) = sum(Y)
%         sums(2) = sum(Y.^2)
%         where Y are iid samples with expected value:
%         E[P_0]           on level 0
%         E[P_l - P_{l-1}] on level l>0
%         cost = cost of N samples

r=1;

sums = zeros(1, 7);
pi_f = zeros(1, N);
pi_c = zeros(1, N); 
for i = 1:N

    if l == 0
            nf = 8.^l;
            xf = 2*rand(1,nf)-1;      % Fine samples between -1 and 1
            yf = 2*rand(1,nf)-1;
            r2f = (xf).^2 + (yf).^2; % Fine samples radius
            pi_f = sum(r2f<=r^2)*4./(nf); %
            pi_c = 0;
    else
            nf = 8.^l; % Number of fine samples based on level l
            nc = 8.^(l-1);
            xf = 2*rand(1,nf)-1;      % Fine samples between -1 and 1
            yf = 2*rand(1,nf)-1;      % Fine samples between -1 and 1
            
            xc = xf(1:nc);
            yc = yf(1:nc); 
            r2c = xc.^2 + yc.^2;       % compute squared distance to origin
            r2f = (xf).^2 + (yf).^2; 
            
            pi_f = sum(r2f<=r^2)*4./(nf); % the area is 4 (area of square) times the proportion of points in the circle. 
            pi_c = sum(r2c<=r^2)*4./(nc); 
    end

        sums(1) = sums(1) + sum(pi_f-pi_c); % P_f - P_l
        sums(2) = sums(2) + sum((pi_f-pi_c).^2); % Variance
        sums(3) = sums(3) + sum((pi_f-pi_c).^3); % Skewness
        sums(4) = sums(4) + sum((pi_f-pi_c).^4); % Kurtosis 
        sums(5) = sums(5) + sum(pi_f);
        sums(6) = sums(6) + sum(pi_f.^2); 
        sums(7) = sums(7) + sum(pi_c);
        
end

cost = nf*N;
