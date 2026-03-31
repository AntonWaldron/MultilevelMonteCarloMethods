% This function produces the Brownian motion and Brownian Bridging plots used in the FYP. 
function BPATH1

figure
randn('state',100) % set the state of randn
T = 1; N = 500; dt = T/N; t = [dt:dt:1];

M = 1000;  % M paths simultaneously
dW = sqrt(dt)*randn(M,N); % increments
W = cumsum(dW,2); % cumulative sum
Wmean = mean(W); 
plot([0,t],[1,Wmean],'b-'), hold on % plot mean m over M paths


h1 = plot(t, Wmean, 'k', 'LineWidth', 1);          % mean of 1000 paths
h2 = plot([0,t],[ones(10,1), W(1:10,:)],'r--');     % 10 individual paths

xlabel('t','FontSize',16)
ylabel('W(t)','FontSize',16,'Rotation',0,'HorizontalAlignment','right')

legend([h1 h2(1)], ...
       {'mean of 1000 paths','10 brownian motion paths'}, ...
       'Location','northwest')

hold off

figure
% Brownian Bridge Construction
T = 1; 
N = 500; 
dt = T/N; 
t = (dt:dt:T)';  % column vector (N x 1)

M = 1000;
dW = sqrt(dt)*randn(N,M);   % (N x M) increments
W = cumsum(dW,1);           % (N x M) standard BM paths

% Brownian bridge
W_bridge = W - (t) * W(end,:); % (N x M) each path ends at 0
Wmean = mean(W_bridge, 2);
size(Wmean)
size(t)
plot(t, Wmean), hold on
plot(t,W_bridge(:,1:10), 'r--'), hold off
xlabel('t', 'FontSize', 16)
ylabel('W(t)', 'FontSize', 16, 'Rotatio', 0, 'HorizontalAlignment', 'right')
legend('Mean of 1000 paths', '10 Individual Brownian Bridges')
