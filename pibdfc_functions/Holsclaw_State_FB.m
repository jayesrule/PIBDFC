function [states, trans_prob] = Holsclaw_State_FB(sqrl, Z, xDat, Omega, X)
%% Inputs:
% Xi : parameters controlling the transition probabilities and the log linear regression term
% Z : Matrix with the Markov property, binary matrix showing which states are present at each time point
% Exo : Exogenous variables used in the transition probabilities
% Omega : The precision matrices associated with each class
% Y : The observations

%% Outputs
% states : sampled states for current iteration
% marg_prob : marginal probability of states at each time point
% trans_prob : transition probabilities

n = size(xDat,1);
H = size(xDat,2);
k = size(Z,2);
p = size(Omega,1);
% Make Exo the concatenation of Z and X, so that Exo \in \mathbb{R}^{T\times K+B} where $B$ is the number
% of covariates
Exo = cat(2,Z,xDat)';
mu = repelem(0,p);

% Storage
trans_prob = zeros(k,k,n);
states = zeros(1,n);
probs = zeros(1,k);
% Need to find a way to vectorize this
%%
for t = 1:n 
    for i = 1:k 
        trans_prob(i,:,t) = exp(sqrl(i,(1:k))+(sqrl(:,k+(1:H))*Exo(k+(1:H),t))');
    end
    trans_prob(:,:,t) = trans_prob(:,:,t)./repmat(sum(trans_prob(:,:,t),2),1,k);
end

%Initializing P for forward Recursion
p = zeros(k,k,n);

%Keeping track of log pi in the recursion
lpi = zeros(n,k);

%Getting pi1
l = zeros(k,1);
x = X(1,:);
for i = 1:k
    l(i) = log(trans_prob(i,find(Z(2,:)),2))+....
                mvnloglike(X(1,:),mu,Omega(:,:,i));
end

lstar = sum(exp(l));
lpi(1,:) = l-log(lstar);

%%
%Begin Forward Recursion
for v = 2:n
    for r = 1:k
        for s = 1:k
            %log-liklihood of being in state s
            loglik = mvnloglike(X(v,:),mu,Omega(:,:,s));
            p(r,s,v) = lpi(v-1,r)+log(trans_prob(r,s,v))+loglik;
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
states = repelem(0,n);

%Sampling state at the nth time point
states(n) = samp_state(pis(n,:));

%Backward sampling of states
for i = 1:(n-1)
    states(n-i) = samp_state(p(:,states(n-i+1),n-i+1)./sum(sum(p(:,states(n-i+1),n-i+1))));
end

