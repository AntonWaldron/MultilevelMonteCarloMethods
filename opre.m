%
% These are similar to the MLMC tests for the original 
% 2008 Operations Research paper, using an Euler-Maruyama 
% discretisation with 4^l timesteps on level l.
%
% The differences are:
% -- the plots do not have the extrapolation results
% -- the top two plots are log_2 rather than log_4
% -- the new MLMC driver is a little different
% -- switch to X_0=100 instead of X_0=1
%

function opre

close all; clear all;

addpath('..');

N0    = 1000;   % initial samples on coarse levels
Lmin  = 2;      % minimum refinement level
Lmax  = 6;      % maximum refinement level
for option = 1:6
  if (option==1) 
    fprintf(1,'\n ---- European call ---- \n');
    N      = 1000;       % samples for convergence tests
    L      = 5;             % levels for convergence tests 
    Eps    = [ 0.005 0.01 0.02 0.05 0.1 ];
  elseif (option==2) 
    fprintf(1,'\n ---- Asian call ---- \n');
    N      = 1000;  % samples for convergence tests
    L      = 5;        % levels for convergence tests 
    Eps    = [ 0.005 0.01 0.02 0.05 0.1 ];
  elseif (option==3) 
    fprintf(1,'\n ---- lookback call ---- \n');
    N      = 1000;  % samples for convergence tests
    L      = 5;        % levels for convergence tests 
    Eps    = [ 0.01 0.02 0.05 0.1 0.2 ];
  elseif (option==4) 
    fprintf(1,'\n ---- digital call ---- \n');
    N      = 4000;  % samples for convergence tests
    L      = 5;        % levels for convergence tests 
    Eps    = [ 0.02 0.05 0.1 0.2 0.5 ];
  elseif (option==5) 
    fprintf(1, '\n ---- Fixed Lookback call ---- \n')
      N      = 2000;  % samples for convergence tests
      L      = 5;        % levels for convergence tests 
      Eps    = [ 0.02 0.05 0.1 ];
  elseif (option==6)
      fprintf(1,'\n ---- Heston model ---- \n');
        N      = 2000;  % samples for convergence tests
        L      = 5;        % levels for convergence tests 
        Eps    = [ 0.01 0.02 0.05 0.1 ];
  end

  if (option<6)
    filename = ['option' num2str(option)];
  else
    filename = 'opre_heston_nobs';
  end
  fp = fopen([filename '.txt'],'w');
  mlmc_test_2(@opre_l, N,L, N0,Eps,Lmin,Lmax, fp, option);
  fclose(fp);

%
% plot results
% % Changed to pdf from eps
  nvert = 3;
  mlmc_plot_2(filename, nvert);

  if(nvert==1)
    figure(1)
    print('-dpdf',[filename 'a.pdf'])
    figure(2)
    print('-dpdf',[filename 'b.pdf'])
  else
    print('-dpdf',[filename '.pdf'])
  end

%
% print exact analytic value, based on S0=K
%
  T   = 1;
  r   = 0.05;
  sig = 0.2;
  K   = 100;
  S0 = 100;
  k   = 0.5*sig^2/r;
  d1  = (r+0.5*sig^2)*T / (sig*sqrt(T));
  d2  = (r-0.5*sig^2)*T / (sig*sqrt(T));
  
  % Fixed Lookback pricing formula addition - AW
  b1 = (-r + sig.^2./2)/(sig.*1); 
  b2 = b1 - sig;
  b3 = (r - sig.^2./2)./sig;
  Y2 = 0;
  p = S0*exp(-r*T)*(normcdf(b1) - (sig.^2./(2.*r)).*exp(Y2)*normcdf(-b3)) + S0.*(sig.^2./(2.*r)).*normcdf(-b2) - S0.*normcdf(b2);


  if (option==1)
    val = NaN;
  elseif (option==2)
    val = NaN;
  elseif (option==3)
    val = K*( ncf(d1) - ncf(-d1)*k - exp(-r*T)*(ncf(d2) - ncf(d2)*k) );
  elseif (option==4)
    val = K*exp(-r*T)*ncf(d2);
  elseif (option==5)
    val = p + S0 - K.*exp(-r.*T);
  elseif (option==6)
      val=NaN;
  end

  if isnan(val)
    fprintf(1,'\n Exact value unknown \n\n');
  else
    fprintf(1,'\n Exact value: %f \n\n',val);
  end

%
% now do 100 MLMC calcs in parallel
%
  filename = [filename 'fixed_n'];
  fp = fopen([filename '.txt'],'w');
  % Had to change from mlmc_test_n to 100
  mlmc_test_100(@opre_l, val, N0,Eps,Lmin,Lmax, fp, option);
  fclose(fp);

%
% plot results
% had to change from mlmc_plot_n to 100
  mlmc_plot_100(filename);
  print('-dpdf',[filename '.pdf'])
end
end



%-------------------------------------------------------
%
% level l estimator for Operations Research paper
%

function [sums, cost] = opre_l(l,N, option)

M = 4; % The factor by which the timestep is refined at each level

T   = 1;
r   = 0.05;
sig = 0.2;
K   = 100;

nf = M^l;    % Number of timesteps on finest level
nc = M^(l-1); %nf/M;   % Number of timesteps on coarse level
 
hf = T/nf;   % Timestep for the fine level
hc = T/nc;   % Timestep for the coarse level


sums(1:6) = 0;

for N1 = 1:10000:N 
  N2 = min(10000,N-N1+1);  % This is always 10000 ? 

%
% Geometric Brownian Motion model i.e. not Heston Model
%
  if option<6

    % Set initial price to strike price  
    X0 = K;
    
    % Initialise coarse and fine to represent all paths starting at X_0 =
    % K
    Xf = X0*ones(1,N2);
    Xc = Xf;
    
    % Asian Options are estimated using the trapezoidal rule, so we start with
    % Half weighted initial value 1/2X_{0}h
    Af = 0.5*hf*Xf;
    Ac = 0.5*hc*Xc;

    % For Lookback Option
    Mf = Xf;
    Mc = Xc;

    % For Fixed Lookback max
    mf = Xf;
    mc = Xc;
    
    % For level 0 
    if l==0
         dWf = sqrt(hf)*randn(1,N2); % A single fine-level Brownian increment
         Xf  = Xf + r*Xf*hf + sig*Xf.*dWf; % Initial estimate of stock price using Euler-Marayama
         Af = Af + 0.5*hf*Xf; % First half of the trapezoidal rule for Asian Option
         Mf = min(Mf,Xf); % Initial price minimum for Lookback Option
         mf = max(mf, Xf);
    else
        % For levels > 0 nc timesteps
      for n = 1:nc
        % Initialising the Brownian increment for coarse timesteps 
        dWc = zeros(1,N2);
        for m = 1:M % Factor by which the timestep is refined at each level
          
          % Brownian motion term for the fine levels
          dWf = sqrt(hf)*randn(1,N2);

          % Brownian motion term for the coarse levels are dependent on the
          % fine terms to ensure coupling 
          dWc = dWc + dWf; 
            
          % Euler Marayama for approximating SDE solution X. 
          Xf  = Xf + r*Xf*hf + sig*Xf.*dWf;
          % Asian Option: Updating the average of the stock price 
          Af  = Af + hf*Xf; % hfXf = Xf/nf 
          
          % Lookback Option: Updating the minimum of the stock price
          Mf  = min(Mf,Xf);
          mf = max(mf, Xf); % Updating max
        end

        % Euler Marayama for coarse terms using dWc = dWc + dWf to ensure coupling  
        Xc = Xc + r*Xc*hc + sig*Xc.*dWc;
        % Average stock price over the coarse estimates
        Ac = Ac + hc*Xc;

        % Minimum of the stock price, over the coarse terms for Lookback
        Mc = min(Mc,Xc);

        mc = max(mc, Xc); 
      end

      % Final X_{f}, X_{c} values had full weight rather than half. 
      Af = Af - 0.5*hf*Xf;
      Ac = Ac - 0.5*hc*Xc;
    end

    % Various option payoffs are now computed
    if option==1
      Pf = max(0,Xf-K);
      Pc = max(0,Xc-K);
    elseif option==2
      Pf = max(0,Af-K);
      Pc = max(0,Ac-K);
    elseif option==3 % floating lookbacks
      beta = 0.5826;  % special factor for offset correction
      Pf = Xf - Mf*(1-beta*sig*sqrt(hf));
      Pc = Xc - Mc*(1-beta*sig*sqrt(hc));
    elseif option==4 % Digital option
      Pf = K * 0.5*(sign(Xf-K)+1);
      Pc = K * 0.5*(sign(Xc-K)+1);
    elseif option == 5 % Adding Fixed Lookback call option
      beta = 0.5826;
      Pf = max(0, mf*(1+beta*sig*sqrt(hf)) - K);
      Pc = max(0, mc*(1+beta*sig*sqrt(hc)) - K);
    end
    % 

    % This else loop relates to the heston model
  else

    if option == 7
        X0 = [K; 0.04];
        Xf = X0*ones(1,N2);
        Xc = Xf;
    
        if l==0
          dWf = sqrt(hf)*randn(2,N2);
          Xf  = Xf + mu(Xf,hf)*hf + sig_dW(Xf,dWf,hf);
    
        else
          for n = 1:nc
            dWc = zeros(2,N2);
            for m = 1:M
              dWf = sqrt(hf)*randn(2,N2);
              dWc = dWc + dWf;
              Xf  = Xf + mu(Xf,hf)*hf + sig_dW(Xf,dWf,hf);
            end
            Xc = Xc + mu(Xc,hc)*hc + sig_dW(Xc,dWc,hc);
          end
        end
    
        Pf = max(0,Xf(1,:)-K);
        Pc = max(0,Xc(1,:)-K);
    end
   
  end
    % Heston model ends
    
  % Discounting the expected option cost   
  Pf = exp(-r*T)*Pf;
  Pc = exp(-r*T)*Pc;

  if l==0
    Pc=0;
  end

  sums(1) = sums(1) + sum(Pf-Pc); % P Sure [P_l - P_l-1]

  % A sum of the above is used for the total MLMC estimator.

  % Then below are used for like kurtosis and variance I think: 
  sums(2) = sums(2) + sum((Pf-Pc).^2); % sums(2) = sums(Y.^2) 
  sums(3) = sums(3) + sum((Pf-Pc).^3); % ..etc
  sums(4) = sums(4) + sum((Pf-Pc).^4);
  sums(5) = sums(5) + sum(Pf); 
  sums(6) = sums(6) + sum(Pf.^2); 
end

cost = N*nf;   % cost defined as number of fine timesteps * Number of samples
end

%--------------------

function m=mu(x,h)

%m = [ 0.05*x(1,:); ...
%       5*(0.04-x(2,:)) ];

m = [ 0.05*x(1,:); ...
       ((1-exp(-5*h))/h)*(0.04-x(2,:)) ];
end

%--------------------

function sigdW=sig_dW(x,dW,h)

dW(2,:) = -0.5*dW(1,:) + sqrt(0.75)*dW(2,:);

%sigdW = [ sqrt(max(0,x(2,:))).*x(1,:).*dW(1,:);  ...
%          0.25*sqrt(max(0,x(2,:))).*dW(2,:) ];

sigdW = [ sqrt(max(0,x(2,:))).*x(1,:).*dW(1,:);  ...
          exp(-5*h)*0.25*sqrt(max(0,x(2,:))).*dW(2,:) ];
end

%
% Normal CDF function
%

function N = ncf(x)

N = 0.5*erfc(-x/sqrt(2));

end
