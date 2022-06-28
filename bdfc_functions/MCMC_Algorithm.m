function [C_save, Sig_save, adj_save, Theta_save, ar_gamma, ar_theta, ...
    nu_save, trans_save, ar_nu, ppi_HMM, states,edge_probs,...
	states_ar,states_prop_storage, scott_ar, base_ar,num_prop_scott, num_prop_base]...
    = MCMC_Algorithm(V0, V1, lambda, pii, alpha, beta,my_w, sigma_alpha_1,...
	alpha_0, burnin, nmc, disp, yDat, xDat, S)
	
% Modified from code provided by Christine Peterson and Hao Wang to allow inference on multiple graphs
% Input parameters:
%   Theta : initial value for graph similarity matrix
%   b_prior : vector of degrees of freedom parameters for G-Wishart priors for each group
%   D_prior : p x p x K array of location parameters for G-Wishart prior for each group
%   n: vector of sample sizes for each group
%   S : p x p x K array of sample covariance matrix for each group
%   C : p x p x K array of initial precision matrices for each group
%   nu : p x p matrix of initial values for nu
%   alpha : shape parameter for Theta slab
%   beta : rate parameter for Theta slab
%   a : first parameter for prior on nu
%   b : second parameter for prior on nu
%   my_w : Bernoulli prior parameter for graph similarity indicators
%   burnin : number of burn-in iterations
%   nmc : number of MCMC sample iterations saved
%   disp : T/F for whether to print progress to screen
%   yDat : matrix containing the response
%   xDat : matrix containing the covariates
%   Sdwell : vector containing the length of each block
% Output parameters:
%   C_save : p x p x K x nmc sample of precision matrix
%   Sig_save : p x p x K x nmc sample of covariance matrix
%   adj_save : p x p x K x nmc sample of adjacency matrix
%   Theta_save : K x K x nmc sample of graph similarity matrix
%   ar_gamma : acceptance rate for between-model moves
%   ar_theta : acceptance rate for within-model moves
%   nu_save : p x p x nmc sample of edge-specific parameters
%   ar_nu : acceptance rate for nu
%   lambda_post : posterior samples of lambda
%   beta_post : posterior samples of beta
%   accept_lambda : acceptance rate of lambda
%   accept_beta : acceptance rate of beta
%   ppi_HMM : posterior probability of inclusion for the HMM classes
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=size(xDat,1);
p = size(yDat, 2);
T=size(xDat,2);
%display(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Proposal parameters for MH steps of the graph procedure
alpha_prop = 1;
beta_prop = 1;

tau = V1;

ind_noi_all = zeros(p-1,p);
 
for i = 1:p
       if i==1  
       ind_noi = [2:p]'; 
      elseif i==p
       ind_noi = [1:p-1]'; 
      else
       ind_noi = [1:i-1,i+1:p]';
       end
       ind_noi_all(:,i) = ind_noi;
end

pii_RB = zeros(p, p, S);
pii_mat = zeros(p, p, S);
for i = 1:S
    pii_mat(:,:,i) = pii(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Posterior Storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_save = zeros(p, p, S, nmc);
Sig_save = C_save;
adj_save = C_save;
Theta_save = zeros(S, S, nmc);
nu_save = zeros(p, p, nmc);
mu_save=zeros(p,nmc);
trans_save = zeros(S,S,T,nmc);
edge_probs = zeros(p,p,S,nmc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceptance Rate Stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ar_gamma = zeros(S, S);
n_within_model = zeros(S, S);
ar_theta = zeros(S, S);
ar_nu = zeros(p, p);
accept_lambda=zeros(1,p); %acceptance count for lambda
accept_gamma_beta=zeros(p,K); %Acceptance for variable selection with respect to beta
ppi_HMM=zeros(S,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adj = repmat(eye(p),1,1,S);
Sig=repmat(eye(p),1,1,S);
C = Sig;
%for s = 1:S 
%   Sig(:,:,s)=cov(yDat(logical(xDat(s,:)),:),1)';
%end

Theta = zeros(S,S);
nu = zeros(p,p);
for i = 1:p
    for j = 1:p
        nu(i,j) = alpha_0;
    end
end
Y=yDat';
states_Init = randsample(1:S,T,true);
store = zeros(S,T);
for t = 1:T
    store(states_Init(t),t)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HIDDEN MARKOV MODEL STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
states=zeros(T,nmc);
%Initialize the matrix with the state counts for HMM
st_counts=zeros(S,S);

% Define logit function for Sharon prior on structural T1 data
logit = @(x) exp(x)./(1+exp(x));

% Initialize alpha_1
alpha_1 = 0;
ar_alpha = 0;

% Acceptance rate stuff for the HMM sampler
states_ar = 0;
num_prop_scott= 0;
scott_ar = 0;
num_prop_base = 0;
base_ar = 0;

% Storage for HMM sampler objects
fp_store = zeros(S,S,T);
bp_store = zeros(S,S,T);

% Storage to return all proposed states
states_prop_storage = zeros(T,burnin+nmc);

%Create initial array of covariance matrices
Scov=Y(:,store(1,:)==1)*Y(:,store(1,:)==1)';
for s = 2:S    
	Scov=cat(3,Scov,Y(:,store(s,:)==1)*Y(:,store(s,:)==1)');
end


%%%%%% THIS ALL NEEDS TO BE EXPORTED TO THE CALL FUNCTION %%%%%%
a = 2;
b = 2;
c = 0;
d = 10;
a_prop = 2;
b_prop = 2;

Xi = zeros(S,S+size(xDat,1)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This step is to compute Sig_baseline, only performed once
%display(size(store))
n = sum(store,2)';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform MCMC sampling
for iter = 1:burnin + nmc
    if (disp && mod(iter, 1000) == 0)
        fprintf('iter = %d\n', iter);
        sum(store',1)
    end

%% MH and Gibbs sampling step for Hidden Markov Model
   % Update graph and precision matrix for each group using code from Hao Wang
    for cur_graph = 1:S

		% Sample Sig and C
		for i = 1:p
			ind_noi = ind_noi_all(:,i);


			tau_temp = tau(ind_noi,i,cur_graph);

			Sig11 = Sig(ind_noi,ind_noi,cur_graph); 
			Sig12 = Sig(ind_noi,i, cur_graph);

			invC11 = Sig11 - Sig12*Sig12'/Sig(i,i, cur_graph);

			Ci = (Scov(i,i,cur_graph)+lambda)*invC11+diag(1./tau_temp);  


			Ci = (Ci+Ci')./2;       
			Ci_chol = chol(Ci);    
			mu_i = -Ci_chol\(Ci_chol'\Scov(ind_noi,i,cur_graph));
			epsilon = mu_i+ Ci_chol\randn(p-1,1);


			C(ind_noi, i, cur_graph) = epsilon;
			C(i, ind_noi, cur_graph) = epsilon;

			a_gam = 0.5*n(:,cur_graph)+1;
			b_gam = (Scov(i,i,cur_graph)+lambda)*0.5;
			gam = gamrnd(a_gam,1/b_gam);

			c = epsilon'*invC11*epsilon;

			C(i,i, cur_graph) = gam+c;

			% note epsilon is beta in original Hao Wang code


			%% Below updating Covariance matrix according to one-column change of precision matrix
			invC11epsilon = invC11*epsilon;

			Sig(ind_noi,ind_noi,cur_graph) = invC11+invC11epsilon*invC11epsilon'/gam;
			Sig12 = -invC11epsilon/gam;
			Sig(ind_noi,i,cur_graph) = Sig12;
			Sig(i,ind_noi,cur_graph) = Sig12';
			Sig(i,i,cur_graph) = 1/gam;


			v0 = V0(ind_noi,i,cur_graph);
			v1 = V1(ind_noi,i,cur_graph);

            mrf_sum=0;
            for foo = 1:S
                mrf_sum=mrf_sum+Theta(cur_graph,foo)*pii_mat(ind_noi,i,foo);
            end

            % Applying the MRF prior to probability of inclusion
            pii_mat(ind_noi,i,cur_graph)=exp(pii_mat(ind_noi,i,cur_graph).*(nu(ind_noi,i)+2*mrf_sum))./(1+exp(nu(ind_noi,i)+2*mrf_sum));     

            w1 = -0.5*log(v0) -0.5*epsilon.^2./v0+log(1-pii_mat(ind_noi,i,cur_graph));
            w2 = -0.5*log(v1) -0.5*epsilon.^2./v1+log(pii_mat(ind_noi,i,cur_graph));

            w_max = max([w1,w2],[],2);
            w = exp(w2-w_max)./sum(exp([w1,w2]-repmat(w_max,1,2)),2);
			 
            pii_mat(ind_noi,i,cur_graph) = w;
            pii_mat(i, ind_noi,cur_graph) = w;
            z = (rand(p-1,1)<w);

            v = v0;
            v(z) = v1(z);

            tau(ind_noi, i, cur_graph) = v;        
            tau(i, ind_noi, cur_graph) = v;

            adj(ind_noi, i, cur_graph) = z;
            adj(i, ind_noi, cur_graph) = z;

		end
		
        if iter > burnin
			edge_probs(:,:,cur_graph,iter-burnin) = pii_mat(:,:,cur_graph);
            C_save(:, :, cur_graph, iter-burnin) = C(:, :, cur_graph);
            adj_save(:, :, cur_graph, iter-burnin) = adj(:, :, cur_graph);
            Sig_save(:,:,cur_graph, iter-burnin) = Sig(:,:,cur_graph); 

        end
    end
%% Holsclaw et al. HMM Sampler
    % Update parameters
    Xi = Holsclaw_Param_Vec(Xi, store', xDat(2:end,:)', c, d);
    % Update states and probabilities
    [states_Init, trans_prob] = Holsclaw_State2(Xi, store', xDat(2:end,:), Sig, yDat);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Given sampled states create binary indicator mtarix of state classification
    for s = 1:S
        for t = 1:T 
            if(states_Init(t)==s);
                store(s,t)=1;
				% If burn in reached add to posterior probability for each 
                if iter>burnin
                   ppi_HMM(s,t)=ppi_HMM(s,t)+1; 
                end
            else
                store(s,t)=0;
            end
        end
    end 

	% If burn in reached then add to state posterior storage
    if iter>burnin
        states(:,iter-burnin)=states_Init;
    end

    %Obtain matrices for compute the covariance matrices 
    yNorm = Y; 

    %Create array of covariance matrices
    Scov=yNorm(:,store(1,:)==1)*yNorm(:,store(1,:)==1)';
    for s = 2:S    
        Scov=cat(3,Scov,yNorm(:,store(s,:)==1)*yNorm(:,store(s,:)==1)');
    end
	
	n = sum(store,2)';
		
% 	if mod(iter,10)==0
% 		display(n)
% 	end
 

% Update the parameters for network relatedness
    for k = 1:S-1
        for m = k+1:S
            % Between model move
            if Theta(k, m) == 0
                theta_prop = gamrnd(alpha_prop, beta_prop);
            else
                theta_prop = 0;
            end
            
            Theta_prop = Theta;
            Theta_prop(k, m) = theta_prop;
            Theta_prop(m, k) = theta_prop;
            
            % Get terms that are a sum over all edges on log scale
            sum_over_edges = 0;
            for i = 1:p-1
                for j = i+1:p
                    sum_over_edges = sum_over_edges + ...
                        log(calc_mrf_C_original(Theta, nu(i, j))) + 2 * ...
                        (theta_prop - Theta(m, k)) * adj(i, j, k) * adj(i, j, m) - ...
                        log(calc_mrf_C_original(Theta_prop, nu(i, j)));
                end
            end
            
            % Calculate MH ratio on log scale
            if theta_prop == 0
                log_ar = alpha_prop * log(beta_prop) - log(gamma(alpha_prop)) + ...
                    log(gamma(alpha)) - alpha * log(beta) - ...
                    (alpha - alpha_prop) * log(Theta(m, k)) + ...
                    (beta - beta_prop) * (Theta(m, k)) + sum_over_edges + ...
                    log(1 - my_w) - log(my_w);
            else
                log_ar = alpha * log(beta) - log(gamma(alpha)) + ...
                    log(gamma(alpha_prop)) - alpha_prop * log(beta_prop) - ...
                    (alpha - alpha_prop) * log(theta_prop) - ...
                    (beta - beta_prop) * theta_prop + sum_over_edges + ...
                    log(my_w) - log(1 - my_w);
            end
            
            % Accept proposal with given probability
            if log_ar > log(unifrnd(0,1))
                Theta(m, k) = theta_prop;
                Theta(k, m) = theta_prop;
                
                % Increment acceptance rate
                ar_gamma(k, m) = ar_gamma(k, m) + 1 / (burnin + nmc);
            end
            
            % Within model model
            if Theta(k, m) ~= 0
                n_within_model(k, m) = n_within_model(k, m) + 1;
                theta_prop = gamrnd(alpha_prop, beta_prop);
                Theta_prop = Theta;
                Theta_prop(k, m) = theta_prop;
                Theta_prop(m, k) = theta_prop;
                
                % Get terms that are a sum over all edges on log scale
                sum_over_edges = 0;
                for i = 1:p-1
                    for j = i+1:p
                        sum_over_edges = sum_over_edges + ...
                            log(calc_mrf_C_original(Theta, nu(i, j))) + 2 * ...
                            (theta_prop - Theta(m, k)) * adj(i, j, k) * adj(i, j, m) - ...
                            log(calc_mrf_C_original(Theta_prop, nu(i, j)));
                    end
                end
                
                % Calculate MH ratio on log scale
                log_theta_ar = (alpha - alpha_prop) * (log(theta_prop) - log(Theta(m, k))) + ...
                    (beta - beta_prop) * (Theta(m, k) - theta_prop) + sum_over_edges;
                
                % Accept proposal with given probability
                if log_theta_ar > log(unifrnd(0,1))
                    Theta(m, k) = theta_prop;
                    Theta(k, m) = theta_prop;
                    
                    % Track number of proposals accepted
                    ar_theta(k, m) = ar_theta(k, m) + 1;
                end
            end
        end
	end
	
%%  Update Marginal Edge Parameters 
    % Generate independent proposals for q from beta(a_prop, b_prop) density
    for i = 1:p-1
        for j = i+1:p
            q = betarnd(a_prop, b_prop);
            nu_prop = log(q) - log(1-q);
            
            % Calculate MH ratio on log scale
            % log(p(nu_prop)) - log(p(nu)) + log(q(nu)) - log(q(nu_prop))
            log_nu_ar = (nu_prop - nu(i, j)) * (sum(adj(i, j, :)) + a - a_prop) - ...
                (a + b - a_prop - b_prop) * log(1 + exp(nu_prop)) - ...
                log(calc_mrf_C_original(Theta, nu_prop)) + ... 
                (a + b - a_prop - b_prop) * log(1 + exp(nu(i, j))) + ...
                log(calc_mrf_C_original(Theta, nu(i, j)));
            
            if log_nu_ar > log(unifrnd(0,1))
                nu(i, j) = nu_prop;
                nu(j, i) = nu_prop;
                ar_nu(i, j) = ar_nu(i, j) + 1 / (burnin + nmc);
            end
        end
    end
    
    % Retain values for posterior sample
    if iter > burnin
        C_save(:, :, :, iter-burnin) = C;
        adj_save(:, :, :, iter-burnin) = adj;
        nu_save(:, :, iter-burnin) = nu;
        Theta_save(:, :, iter-burnin) = Theta;
        Sig_save(:, :, :, iter-burnin) = Sig;
        trans_save(:,:,:,iter-burnin) = trans_prob;
    end
end %END ITERATIONS

%%
% Compute acceptance rate for theta as number of proposals accepted divided
% by number of within model moves proposed
for k = 1:S-1
    for m = k+1:S
        ar_theta(k, m) = ar_theta(k, m) / n_within_model(k, m);
    end
end
edge_probs = sum(edge_probs.*adj_save,4)./sum(adj_save,4);

for s = 1:S
	edge_probs(:,:,s) = edge_probs(:,:,s).*(~eye(p)) + eye(p);
end
states_ar = states_ar / nmc;
base_ar = base_ar/num_prop_base;
scott_ar = scott_ar / num_prop_scott;
%%
% Normalize acceptance ratios and ppi for classification
accept_lambda = accept_lambda/(nmc);
accept_gamma_beta=accept_gamma_beta/(nmc);
ppi_HMM=ppi_HMM/(nmc);


