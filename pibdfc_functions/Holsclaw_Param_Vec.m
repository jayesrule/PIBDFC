function [sqrl] = Holsclaw_Param_Vec(sqrl_old, st_ind, xDat, a, b)

%% This function is a Gibbs sampling step to sample new values for the log linear coefficients 
%% in the Holsclaw et al. paper for sampling from heterogenous HMMs.

%% Inputs:
% Xi_old: The values of the Xi parameter for the previous iteration
% Z: The arrangement of the state transitions in binary form for the current iteration
% Exo: The exogenous variables incorporated in the log linear regression, X in the Holsclaw paper
% a: prior mean of the Xi, can be scalar or assign unique prior value to each term
% b: prior variance of the Xi, can be scalar or assign unique prior value to each term

%% Outputs:
% Xi_new: The sampled values of the Xi parameter for the current iteration

%%%%%%
%% Polya Gamma sampler is out of date and can only sample scalars, the Holsclaw algorithm currently
%% samples each omega conditional on all others, so this isn't an issue, but it's something to fix if 
%% I ever want to use this again.
%% Actually I need to fix this
%%%%%%

% Make Exo the concatenation of Z an X, so that Exo \in \mathbb{R}^{T\times K+B} where $B$ is the number
% of covariates


%%
B = size(st_ind,2);
%Z = cat(2,zeros(B,1),st_ind(2:(end),:)')';
Z = st_ind;
Exo = cat(2,Z,xDat);

% Index variables
H = size(Exo,2);
K = size(sqrl_old,1);
T = size(Exo,1);

% Expand a and b to approach dimension if not given as predictor specific
try
    a = a.*ones(size(sqrl_old));
catch
    warning('a needs to be a matrix of dimension equal to the covariate dimension, or a scalar');
end
try
    b = b.*ones(size(sqrl_old));
catch
    warning('b needs to be a matrix of dimension equal to the covariate dimension, or a scalar');
end
%%
% Storage for variables
Omega_kh = zeros(T,T);
eta_kh = zeros(1,T);
C_kh = zeros(T,1);
V_kh = zeros(K,K);
m_kh = 0;
sqrl = sqrl_old;
for(k = 2:K)
    ind_nok = setdiff(1:K,k);
    % Maybe precompute this so it doesn't get evaluated for every MCMC iteration
    for(h = 1:H)
        for(t = 1:T)
             C_kh(t) = log(sum(exp(Exo(t,h).*sqrl(ind_nok,h))));
             eta_kh(t) = Exo(t,h)*sqrl(k,h) - C_kh(t);
             % Once I get pgrnd vectorized delete these twolines
             %Omega_kh(t,t) = pgrnd(eta_kh(t));
        end
        % Need to vectorize pgrnd for this for speed up
        Omega_kh = diag(pgdraw(eta_kh'));
        
        % Mean and Variance of the conditional distribution on $\Xi_{kh}$
        V_kh = (Exo(:,h)'*Omega_kh*Exo(:,h) + b(k,h)^(-1))^(-1);
        m_kh = V_kh*(Exo(:,h)'*((Z(:,k)-.5)-Omega_kh*C_kh)+b(k,h)^(-1)*a(k,h));
        
        % Sample new $\Xi_{kh}$
        sqrl(k,h) = normrnd(m_kh,sqrt(V_kh));
    end
end