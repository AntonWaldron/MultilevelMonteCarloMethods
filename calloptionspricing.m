% This function produces the dataset used for the Multilevel Monte Carlo cost vs Standard cost plot used in Section 3.5 Empirical Results on a Real Life Options Dataset of the FYP. 

function calloptionspricing

close all; clear all;

options_data = readtable('file.csv');

strikePrice = options_data.Strike; % 1x300 column vector of strikeprices
expPrice = options_data.PriceExp; % of stock price at expiry
startPrice= options_data.PriceDate; % Starting stock price
bestBid = options_data.BestBid;
bestOffer = options_data.BestOffer;
target = (bestBid + bestOffer) ./ 2; % Calculate target prices for options
implVol = options_data.ImpliedVolatility; % Implied volatility can be used as variance value
intRate = options_data.Rate./100; % Interest rate in percentage

values = [bestBid, bestOffer, target];
addpath('..');

N0 = 2000; % initial samples on coarse levels
Lmin = 2;  % minimum refinement level
Lmax = 6; % maximum refinement level

% For convergence tests
N = 3000; 
L = 5; 
Eps = [0.01]; 
T = 1./12; % approx 1 month

results = table();
runNames = 1:150;
for i = 1:150
    filename = 'itm_calls_MLMC_eps0_01';
    fp = fopen([filename '.txt'],'w');

    pNames = strcat("P_Run_", string(1:100));
    costNames = strcat("Cost_Run_", string(1:100));
    stdNames = strcat("Std_Run_", string(1:100));
    
    % Now do 100 MLMC calcs in parallel
    % And for 100 MLMC calcs
    fprintf(fp, '\n== 100 MLMC calculations for i = %d ==\n', i);
    R100 = mlmc_test_100(@call_l2, values(i,:), N, Eps, Lmin, Lmax, fp,...
    intRate(i), implVol(i), strikePrice(i), startPrice(i));
    baseData = table(...
        i, ...
        strikePrice(i), ...
        startPrice(i), ...
        implVol(i), ...
        intRate(i), ...
        bestBid(i), ...
        bestOffer(i), ...
        mean([bestBid(i), bestOffer(i)]), ...
     'VariableNames', {'Index', 'StrikePrice', 'StartPrice', 'ImplVol', 'IntRate', 'BestBid', 'BestAsk', 'AvgBidAsk'} ...
    );
    R100_row = [R100(1,:) R100(2,:) R100(3,:)];
    allNames = [pNames costNames stdNames ];
    
    runTable = array2table(R100_row, 'VariableNames', allNames);

    % Combine the metadata and the 100 run results
    rowTable = [baseData runTable];

    % Append to main results table
    results = [results; rowTable];

    fprintf(fp, '---- End of 100 MLMC calcs for i = %d ----\n\n', i);
    % function mlmc_test_100(mlmc_l, val, N0,Eps,Lmin,Lmax, fp, varargin

end


fclose(fp);
outFile = ['mlmc_itm_calls_option_results_' datestr(now,'yyyy-mm-dd_HHMM') '.csv'];
writetable(results, outFile);
fprintf('Results saved to %s\n', outFile);

end



%% Level l estimator for call option with custom inputs

% level l estimator for a call option with inputs

% M^l - number of timesteps used in the simulation
% nf - number of fine timesteps
  
% Inputs:
% l - level used for simulation
% N - number of timesteps
% T - time
% r - interest rate 
% sig - impliedvolatility
% startPrice - starting price of the stock
% strikePrice - K StrikePrice of the option


function [sums,cost] = call_l2(l, N, r, sig, strikePrice, startPrice)
    
M = 4; % The factor by which the timestep is refined at each level
T = 1./12; % These are monthly calls so the timestep is reduced to 1/12
K = strikePrice;

nf = M^l; % Num steps on finest level
nc = M^(l-1); % Num. steps on the coarsest level

hf = T/nf;
hc = T/nc; % Timesteps for each level

sums(1:6) = 0;

for N1 = 1:10000:N
    N2 = min(10000, N-N1+1);

    % Geometric BrownianMotion Model
    X0 = startPrice;
    Xf = X0*ones(1,N2);
    Xc = Xf;

    if l==0
        dWf = sqrt(hf)*randn(1,N2); % A single fine-level Brownian increment
        Xf = Xf + r*Xf*hf + sig*Xf.*dWf; % Initial estimate of stock price using Euler Maruyama
    else

        for n=1:nc
            dWc = zeros(1,N2);
            for m = 1:M

                dWf = sqrt(hf)*randn(1,N2);
                dWc = dWc + dWf; % To ensure coupling between levels
                Xf = Xf + r*Xf*hf + sig*Xf.*dWf; 
            end
            Xc = Xc + r*Xc*hc + sig*Xc.*dWc; % Update coarse-lev
        end
    end


        % Compute option payoff
        Pf = exp(-r*T)*max(0, Xf-K);
        Pc = exp(-r*T)*max(0, Xc-K);

        if l ==0
            Pc=0;
        end

        sums(1) = sums(1) + sum(Pf-Pc);
        sums(2) = sums(2) + sum((Pf-Pc).^2);
        sums(3) = sums(3) + sum((Pf-Pc).^3);
        sums(4) = sums(4) + sum((Pf-Pc).^4);
        sums(5) = sums(5) + sum(Pf);
        sums(6) = sums(6) + sum(Pf.^2);

end
    cost = N*nf;
end



