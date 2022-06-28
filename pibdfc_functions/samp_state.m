function [state] = samp_state(probs)
probs = abs(probs)/sum(abs(probs));
k = length(probs);

state = randsample(1:k,1,true,probs);
end