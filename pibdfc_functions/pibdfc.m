%This function performs the MCMC for pibdfc and stores the output into the
%structure: post_draws

%% Input [T timepoints, R regions, N subjects)
% Y - TxRxN array of the observed data
% xDat - TxKxN array of the exogenous variables
% S - int number of connnectivity states
% Z_prior - SxS matrix of the prior means for the logistic intercept: Z
% eta_prior - SxK matrix of the prior means for the logistic regression
% coef
% a0 - real >0 prior variance the subject level NHMM parameters
% b0 - real >0 prior variance for the population level NHMM parameters
% tau_0 - real >0 prior scale parameter for the graphical horseshoe
% nsims - int number of posterior draws after burn-in
% nburn - int number of posterior draws discarded before saving
% n_report - int report parameter values every n_report samples
% thin - keep only 1/thin of nsims posterior draws
% parallel - boolean should parallel computing be used? Requires parallel
% toolbox

%% output (nsave = nsims/thin; number of posterior samples saved)
% Omega_post - RxRxSxnsave array of posterior samples of Omega
% states_post - TxNxnsave array of posterior samples of state sequences
% sqrl_post - SxSxNxnsave array of posterior samples subject level NHMM
% intercepts
% rho_post - SxKxNxnsave array of posterior samples subject level NHMM
% coefficients
% Z_post - SxSxnsave array of posterior samples group level NHMM intercepts
% eta_post - SxKxnsave array of posterior samples groupt level NHMM
% coefficients
% kappa_post - RxRxSxnsave array of posterior samples of shrinkage
% parameters for Omega




function [post_draws] = pibdfc(Y, xDat, S, Z_prior, eta_prior, a0, b0, tau_0, nsims, nburn, n_report, thin, parallel)


%%
warning('off','MATLAB:nearlySingularMatrix')
[T,R,N] = size(Y);
K = size(xDat,2);


%Creating Posterior storage
nsave = floor(nsims/thin);
Omega_post = zeros(R,R,S,nsave);
states_post = zeros(T,N, nsave);
sqrl_post = zeros(S,S,N,nsave);
rho_post = zeros(S,K,N, nsave);
Z_post = zeros(S,S,nsave);
eta_post = zeros(S,K,nsave);
kappa_post = Omega_post;

%% Setting initial values

Omega = repmat(eye(R),1,1,S);
Omega_inv = Omega;


%Auxillery Graphical Horseshoe parameters
Lambda_sq = repmat(ones(R),1,1,S);
Nu = Lambda_sq;
tau_sq = repelem(1,S);
Xi = tau_sq;


s_assign = ones(T, N);
for i = 1:N
    s_assign(:,i) = randsample(1:S,T,true);
end

%s_assign = state_vec;

%logisitc parameters
sqrl = zeros(S,S,N);
Z = eye(S,S);
Z(1,1) = 0;

rho = zeros(S,K,N);
eta = zeros(S,K);




%% Begining MCMC
for sims = 1:(nsims+nburn)
    %Partition the time series
    part = {};
    for s = 1:S
        part{s} = zeros(sum(sum(s_assign==s)),R);
        ind_count = 1;
        for i = 1:N
            temp = part{s};
            state_count = sum(s_assign(:,i)==s);
            temp(ind_count:(ind_count+state_count-1),:) = Y(s_assign(:,i)==s,:,i);
            part{s} = temp;
            ind_count = ind_count+state_count;
        end
    end
    
    %% Graph
    %Sample the precision parameters
    for i = 1:S
        state_cov = part{i}'*part{i};
        part_lengths(i) = size(part{i},1)+1;
        if part_lengths(i)>2
            [Omega(:,:,i),Omega_inv(:,:,i),...
                Lambda_sq(:,:,i),tau_sq(i), Xi(i),...
                Nu(:,:,i)] = GHS_hmm(state_cov,Omega(:,:,i),...
                Omega_inv(:,:,i),Lambda_sq(:,:,i),...
                tau_sq(i),Nu(:,:,i),...
                Xi(i),part_lengths(i), tau_0);
        end
    end
    
    %% NHMM
    if parallel
        %Sample the state parameters in parallel
        [sqrl, Z, rho, eta] = multi_nhmm_pg_draw_par(xDat, sqrl, Z, Z_prior, rho, eta, eta_prior, s_assign, a0, b0);
        
        %Sample the states
        s_assign = multiHolsclaw_State_draw_par(Y, xDat, Omega, rho, sqrl );
    else
        %sample the state parameters in parallel
        [sqrl, Z, rho, eta] = multi_nhmm_pg_draw(xDat, sqrl, Z, Z_prior, rho, eta, eta_prior, s_assign, a0, b0);
        
        %Sample the states
        s_assign = multiHolsclaw_State_draw(Y, xDat, Omega, rho, sqrl );
        
    end
    
    %% Storage
    if sims > nburn
        % relabel draws
        map = relabel_draws_multi(Omega_prev,s_assign_prev, Omega,s_assign, 0.7);
        
        
        Omega = Omega(:,:,map);
        Omega_inv = Omega_inv(:,:,map);
        Lambda_sq = Lambda_sq(:,:,map);
        tau_sq = tau_sq(map);
        Xi = Xi(map);
        Nu = Nu(:,:,map);
        sqrl = sqrl(map,map,:);
        rho = rho(map,:,:);
        Z = Z(map,map);
        eta = eta(map,:);
        
        lut = uint8([0; map; zeros(255 - size(map, 1), 1)]);
        s_assign = intlut(uint8(s_assign), lut) ;
        
        
        % storing in the output variables
        if (mod(sims,thin) == 0)
            store = floor((sims-nburn)/thin);
            Omega_post(:,:,:,store) = Omega;
            states_post(:,:,store) = s_assign;
            sqrl_post(:,:,:,store) = sqrl;
            rho_post(:,:,:,store) = rho;
            Z_post(:,:,store) = Z;
            eta_post(:,:,store) = eta;
            var_terms = zeros(R,R,S);
            %Calculating the entrywise shrinkage parameter
            for i = 1:S
                var_terms(:,:,i) = tau_sq(i).*Lambda_sq(:,:,i);
            end
            %Storing of shrinkage parameter
            kappa_post(:,:,:,store) = 1./(1+var_terms);
        end
    end
    
    Omega_prev = Omega;
    s_assign_prev = s_assign;
    
    
    %calculate state counts for each subject
    tbl = zeros(N,S);
    for i = 1:N
        tbl(i,:) = histcounts(s_assign(:,i),S);
        if sims < nburn
            for j = 1:S
                if tbl(i,j) < 1.1*R
                    inds = randsample(T,floor(1.1*R));
                    s_assign(inds) = j;
                end
            end
        end
        
    end
    
    
    
    %Progress report
    if mod(sims,n_report)==0
        fprintf("multi nhmm: The iteration number is " + sims)
        tbl
        fprintf('Current parameter positions: \n Z \n')
        Z
        fprintf('Eta \n')
        eta
    end
    
end


%Storing in output variable
post_draws.Omega_post = Omega_post;
post_draws.states_post = states_post;
post_draws.sqrl_post = sqrl_post;
post_draws.rho_post = rho_post;
post_draws.Z_post = Z_post;
post_draws.eta_post = eta_post;
post_draws.kappa_post = kappa_post;
warning('on','MATLAB:nearlySingularMatrix')
end

