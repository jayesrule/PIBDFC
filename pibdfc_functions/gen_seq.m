%This is a function to generate a random sequence of states

function[state_seq] = gen_seq(N, num_states, avg_run)
%%
%%%%%Inputs:
%N (integer): Length of the time series to be returned
%num_states (integer): Number of total states to be considered
%avg_run (integer) : Average run of consecutives state selection.

%%%%Outputs:
%state_seq 1xN time series of states.
%cpoints: vector of change points
%%

states = 1:num_states;

state_seq = double.empty;

i = 1;
n = double.empty;
r = 0;

while r < N
    n = [n,poissrnd(avg_run)];
    i = i+1;
    r = sum(n);
end

n(size(n,2)) = N - sum(n(1:(size(n,2)-1)));
s = randsample(num_states,size(n,2),true);

for j = 1:size(n,2)
    state_seq = [state_seq, repelem(s(j),n(j))];
end

end


