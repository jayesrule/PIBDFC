function [states, trans_prob] = Holsclaw_State(Xi, Z, Exo, Sig, Y)
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

T = size(Exo,2);
H = size(Exo,1);
S = size(Z,2);
P = size(Sig,1);
% Make Exo the concatenation of Z and X, so that Exo \in \mathbb{R}^{T\times K+B} where $B$ is the number
% of covariates
Exo = cat(1,Z',Exo);

% Storage
trans_prob = zeros(S,S,T);
states = zeros(1,T);
probs = zeros(1,S);
% Need to find a way to vectorize this
for(t = 1:T)
    for(i = 1:S)
        trans_prob(i,:,t) = exp(Xi(i,(1:S))+(Xi(:,S+(1:H))*Exo(S+(1:H),t))');
    end
    trans_prob(:,:,t) = trans_prob(:,:,t)./repmat(sum(trans_prob(:,:,t),2),1,S);
end

% Make sure numerical values don't restrict the transition probabilities to be zero and one
%trans_prob = .98*trans_prob + .01;


for(s = 1:S)
    probs(s) = log(trans_prob(s,find(Z(2,:)),2))+....
                logmvnpdf(repmat(Y(1,:),1,1),repmat(0,size(Y(1,:),1),1),Sig(:,:,s));
end
probs = probs - max(probs);
probs = exp(probs)./repmat(sum(exp(probs),2),1,S);
states(1) = find(mnrnd(1,probs',1));

states(1) = 1;
%% Testing with this before implementing the forward backward stuff
for(t = 2:(T-1)) 
    for(s = 1:S)
        probs(s) = log(trans_prob(states(t-1),s,t-1))+...
                    log(trans_prob(s,find(Z(t+1,:)),t+1))+....
                    logmvnpdf(repmat(Y(t,:),1,1),repmat(0,size(Y(t,:),1),1),Sig(:,:,s));
    end
    % Scale probabilities and sample the new state
    probs = probs - max(probs);
    probs = exp(probs)./repmat(sum(exp(probs),2),1,S);
    states(t) = find(mnrnd(1,probs',1));
end

for(s = 1:S)
    probs(s) = log(trans_prob(states(T-1),s,T-1))+...
                logmvnpdf(repmat(Y(T,:),1,1),repmat(0,size(Y(T,:),1),1),Sig(:,:,s));
end

% Scale probabilities and sample the new state
probs = probs - max(probs);
probs = exp(probs)./repmat(sum(exp(probs),2),1,S);
states(T) = find(mnrnd(1,probs',1));
%{
% Stochastic Backward Recursion
for(t = (T-1):1)
    probs = trans_prob(:,states(t+1),t+1);
    states(t) = find(mnrnd(1,probs',1));
end
%}
