%% Function Multilevel Monte Carlo Method for American options

% Inputs: 
%
% l     Level, i..e M^l timesteps used for each path
% N     Number of paths generated 
%
% These two together determine the sizes of the coarse and fine matrices.

function [sums, cost] = america_multilevel(l, N)

T = 1;           % Time to maturity (1 year)
r = 0.06;        % Risk-free rate (6%)
sig = 0.2;       % Volatility (20%)
K = 110;         % Strike price
S0 = 100;        % Initial stock price

M = 4;
if l == 0
    Pc = 0; 
    nf = M^l;% Number of timesteps
    df = T/nf; 
    hf = T/nf;
    dWf = sqrt(df)*randn(nf,N); % N : Number of paths
    Sf  = cumprod([repmat(S0,1,N); exp((r-sig^2/2)*df+sig*dWf)]);  % paths
    Pf = max(K-Sf(end,:),0);  % put option payoff

    % going back in 2 from nf: 
    for n = nf:-1:2
      Pf   = exp(-r*df)*Pf;                            % discount
      itm_f = find(Sf(n,:) < K);  % Filter for in the money paths only
      Xf    = Sf(n,itm_f)';  % Set X and y with itm index for payoff
      Yf    = Pf(itm_f)';
      Af    = [ ones(size(Xf)) (1-Xf) (1-2*Xf+0.5*Xf.^2)]; % 3 basis functions
      beta_f = Af\Yf;                                     % linear regression
      Cf    = Af*beta_f;                                  % continuation value fine                                
      ind_f = itm_f((K-Xf > Cf));                        % find immediate exercise timesteps fine
      Pf(ind_f) = K-Sf(n,ind_f);                            % Payoff due to immediate exercise

    end

else    
    nf = M^l; % Number of timesteps on the fine level
    nc = M^(l-1); % Number of timesteps on the coarse level one below. 
    df = T/nf;
    dc = T/nc;
    hf = T/nf;
    hc = T/nc;
    sums(1:6) = 0; 

    dWf = sqrt(hf)*randn(nf,N);    % (fine timesteps, no. paths)
    dWc = reshape(dWf, M, nc, N);
    dWc = sum(dWc, 1);
    dWc = reshape(dWc, nc, N);
   
    % (coarse timesteps, no. paths)
    Sf  = cumprod([repmat(S0,1,N); exp((r-sig^2/2)*df+sig*dWf)]);  % paths
    Sc  = cumprod([repmat(S0,1,N); exp((r-sig^2/2)*dc+sig*dWc)]);  % paths

    Pf = max(K-Sf(end,:),0);  % put option payoff
    Pc = max(K-Sc(end,:), 0); 
    

    for n = nf:-1:2
      Pf   = exp(-r*df)*Pf;                            % discount
      itm_f = find(Sf(n,:) < K);  % Filter for in the money paths only
      Xf    = Sf(n,itm_f)';  % Set X and y with itm index for payoff
      Yf    = Pf(itm_f)';
      Af    = [ ones(size(Xf)) (1-Xf) (1-2*Xf+0.5*Xf.^2)]; % 3 basis functions
      beta_f = Af\Yf;                                     % linear regression
      Cf    = Af*beta_f;                                  % continuation value fine                                
      ind_f = itm_f((K-Xf > Cf));                % find immediate exercise timesteps fine
      Pf(ind_f) = K-Sf(n,ind_f);  
                                
    end
    
    for n = nc:-1:2
      Pc   = exp(-r*dc)*Pc;                            % discount
      itm_c = find(Sc(n,:) < K);  % Filter for in the money paths only
      Xc    = Sc(n,itm_c)';  % Set X and y with itm index for payoff
      Yc    = Pc(itm_c)';
      Ac    = [ ones(size(Xc)) (1-Xc) (1-2*Xc+0.5*Xc.^2)]; % 3 basis functions
      beta_c = Ac\Yc;                                     % linear regression 
      Cc    = Ac*beta_c;                                  % continuation value fine                                
      ind_c = itm_c((K-Xc > Cc)); 
      Pc(ind_c) = K-Sc(n,ind_c);                            % Payoff due to immediate exercise
    end
   
end


Pf = exp(-r*df)*Pf;

if l == 0
    Pc = 0;
    dc = 0;
end

Pc = exp(-r*dc)*Pc;
    
valf = sum(Pf)/N;
valc = sum(Pc)/N;
   
sums(1) = sum(Pf-Pc); % P Sure [P_l - P_l-1]
sums(2) = sum((Pf-Pc).^2); % sums(2) = sums(Y.^2) 
sums(3) = sum((Pf-Pc).^3); % ..etc
sums(4) = sum((Pf-Pc).^4);
sums(5) = sum(Pf); 
sums(6) = sum(Pf.^2);

cost = N*(nf+nc);
end
