% This is code adapted from 'An Algorithmic Introduction to Numerical Simulation of Stochastic Differential Equations' by Dr. Desmon J. Higham 

% This code produces a plot of how Brownian Motions are coupled in the Multilevel Monte Carlo method. 

function BPATH5

randn('state', 100)
T = 1; N = 256; dt = T/N; t = [dt:dt:1];
M = 1; % 1 path

dWf = sqrt(dt)*randn(M,N); % 1 path x 500 timesteps
dWc = zeros(M, 64); % coarse path

for i = 1:64
        dWc(i) = sum(dWf(4*(i-1)+1:i*4)); 
end

Wf = cumsum(dWf, 2);
Wc = cumsum(dWc,2); 

Nf = length(Wf);
Nc = length(Wc); 

tf = linspace(0, T, Nf);
tc = linspace(0, T, Nc);

plot(tf, Wf', '-x')
hold on
plot(tc, Wc', '-o')
legend('Fine Brownian Motion Path (Wf)', 'Coarse Brownian Motion Path (Wc)')
xlabel('Time')
ylabel('Brownian Motion Value')
title('Brownian Motion Paths')
