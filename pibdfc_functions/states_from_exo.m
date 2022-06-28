function generated_states = states_from_exo(X,k,Q)
%Function to generate states from task data for simulation purposes
%% Input
%X: NxC matrix of task conditions where N is the number of time points and
%C is the number of conditions
%k: the number of states to be drawn
%Q: kxkxC array of transition kernels associated with each condition

%% Output 
%generated_states: Nx1 vector of states drawn from 1:k

%% Declare Variables
[N,C] = size(X);
x = sum(X.*(1:C),2);
generated_states = zeros(N,1);

%Initial state
generated_states(1) = samp_state(repelem(1/k,k));

%Following States
for i = 2:N
    generated_states(i) = samp_state(Q(generated_states(i-1),:,x(i)));
end

end
