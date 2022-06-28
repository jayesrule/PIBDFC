function [states] = multiHolsclaw_State_draw_par(Y, xDat, Omega, rho, sqrl )
%% Inputs:
% Xi : parameters controlling the transition probabilities and the log linear regression term
% Z : Matrix with the Markov property, binary matrix showing which states are present at each time point
% Exo : Exogenous variables used in the transition probabilities
% Omega : The precision matrices associated with each class
% Y : The observations

%% Outputs
% states : sampled states for current iteration

% Storage
k = size(sqrl,1);
[n,R,N] = size(Y);


mu = repelem(0,R);
states = zeros(n,N);
% Need to find a way to vectorize this
%%
parfor sub = 1:N
    
    %Initializing P for forward Recursion
    p = zeros(k,k,n);
    
    %Keeping track of log pi in the recursion
    lpi = zeros(n,k);
    
    %Getting pi1
    l = zeros(k,1);
%     for i = 1:k
%         l(i) = mvnloglike(Y(1,:,sub),mu,Omega(:,:,i));
%     end
%     
%     lstar = sum(exp(l));
%     lpi(1,:) = l-log(lstar);
lpi(1,:) = [0,(zeros(1,k-1)-3)];
    
    %%
    %Begin Forward Recursion
    for v = 2:n
        for r = 1:k
            for s = 1:k
                %log-liklihood of being in state s
                loglik = max(mvnloglike(Y(v,:,sub),mu,Omega(:,:,s)),-400);
                q = sqrl(r,s,sub)+squeeze(xDat(v,:,sub))*squeeze(rho(s,:,sub))'-log(sum(exp(sqrl(r,:,sub)+squeeze(xDat(v,:,sub))*squeeze(rho(:,:,sub))')));
                p(r,s,v) = lpi(v-1,r)+q+0.001+loglik;
            end
        end
        p(:,:,v) = exp(p(:,:,v) - log(sum(sum(exp(p(:,:,v))))));
        lpi(v,:) = log(sum(p(:,:,v),1))';
    end
    
    pis = exp(lpi);
    
    %% Non-stochastic backwards recursion
    PPrime = p;
    piprime = lpi;
    
    for v = flip(2:(n-1))
        for r = 1:k
            for s = 1:k
                piprime(v,:) = log(sum(p(:,:,v+1),1));
                PPrime(r,s,v) = exp(log(p(r,s,v))+piprime(v,s)-lpi(v,s));
            end
        end
    end
    
    
    %% Stochastic Backward Recursion
    temp_states = zeros(n,1);
    
    %Sampling state at the nth time point
    temp_states(n) = samp_state(pis(n,:));
    
    %Backward sampling of states
    for i = 1:(n-1)
        temp_states(n-i) = samp_state(p(:,temp_states(n-i+1),n-i+1)./sum(sum(p(:,temp_states(n-i+1),n-i+1))));
    end
    
    states(:,sub) = temp_states;
    
end