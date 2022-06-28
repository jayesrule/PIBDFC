function [results]...
    =Call_Function(TC, X, S, burnin, nmc, seed)

%% Inputs:
% Y_file :
% Set seed
rng(seed)

%% Call Function Implementation
% Add paths to data
%addpath 'MCMC_Folder'


data.TC = TC;

data.X = X;

% Get the dimensions
p = size(data.TC,2);
T = size(data.TC,1);


base_prob = .5;
alpha_0 =  log(base_prob/(1-base_prob));

% Put in dummy partial correlation for Eric's data without T1 data present
partial_corr = eye(p);

% Slab priors on Theta
param.alpha = 2;
param.beta = 5;
% Bernoulli probability on Theta
param.p_theta = .2;
% Prior variance for alpha_1
param.sigma_alpha_1 = 10;
% Number of states
%S = 4;
% Baseline probability for edge presence in log odds scale
param.alpha_0 = alpha_0;

% Parameters for Hao Wang sampling component of the model
param.h = 50;
param.v0=.02;
param.v1=param.h*param.v0;
param.lambda = 2;
param.pii = 2/(p-1);
param.pii = [param.pii param.pii param.pii param.pii];
param.V0_1 = param.v0*ones(p);
param.V1_1 = param.v1*ones(p);
param.V0 = cat(3,param.V0_1, param.V0_1, param.V0_1, param.V0_1);
param.V1 = cat(3, param.V1_1, param.V1_1, param.V1_1, param.V1_1);
% Correlation matrix (if present) of the T1 structural data for prior
% elicitation
param.structural = partial_corr;

% Get sample sizes for each task from design matrix of stimuli
sampSize=sum(data.X);

% Run the algorithm
% Run the algorithm

%%
[C_save, Sig_save, adj_save, Theta_save, ar_gamma, ar_theta, nu_save, trans_save,ar_nu,...
    ppi_HMM, states, RB_probs,...
    states_ar,states_prop_storage, scott_ar, base_ar,num_prop_scott, num_prop_base] = ...
    MCMC_Algorithm(param.V0, param.V1, param.lambda, param.pii, param.alpha,...
    param.beta, param.p_theta, param.sigma_alpha_1, param.alpha_0, burnin, nmc, true, data.TC, data.X', S);
% PPIs for Theta (graph similarity measure)
ppi_theta = mean(Theta_save ~= 0, 3);

% Edge PPIs for each graph
ppi_edges = mean(adj_save, 4);

% Get 95% credible intervals for omega (precision matrix)
CI_omega_lower = quantile(C_save, 0.025, 4);
CI_omega_upper = quantile(C_save, 0.975, 4);

% Create struct containing the results
w = whos;
for a = 1:length(w)
    results.(w(a).name) = eval(w(a).name);
end

clearvars -except results
