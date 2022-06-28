%This will be the simulation script. Data is generated according to a
%gaussian distribution. Currently requires the Mathematics toolbox to
%generate the precision matrices: Omegas.
addpath('pibdfc_functions');
seed = 2;
rng(seed);


%% Generating the data.
N = 300; %The length of the time series
Subs = 5; %The number of subjects
R = 16; %number of nodes
num_states = 3; %Number of states
%% Creating Block Diagonal Structure with graphs

b_size = [4,2,3]; %size of regions assigned to a block for each state
nblocks = [4,8,6]; %number of blocks in graph. b_size*nblocks >= R

%Initialize Graph
adj = zeros(R,R,num_states);
for i = 1:num_states
    gen = 1*(Generateblkdiag(nblocks(i),b_size(i),1) ~= 0);
    adj(:,:,i) = gen(1:R,1:R);
end


%Setting off-diagonal connections for state 3
adj(10:12,1:3,3) = 1;
adj(1:3,10:12,3) = 1;
adj(7:9,14:16,3) = 1;
adj(14:16,7:9,3) = 1;

%Setting off-diagonal connections for state 2
adj(5:6,9:10,2) = 1;
adj(9:10,5:6,2) = 1;
adj(3:4,7:8,2) = 1;
adj(7:8,3:4,2) = 1;
adj(11:12,15:16,2) = 1;
adj(15:16,11:12,2) = 1;

%%
fprintf('Generating Precision Matrix: Omegas \n')

Omegas = zeros(R,R,num_states); %precision matrix for each state
Sigmas = Omegas; %Covariance matrix for each state

% Generate random covariance
for i = 1:num_states
    eigs = 0;
    
    while(~all(eigs>0)) %Keep generating until S is pos def.
        B = full(sprandsym(adj(:,:,i)));
        O = B'*B-0.5*diag(diag(B)); %Shrink diagonals to boost correlation
        O = 0.5*(O+O'); %symmetrisize
        S = inv(O);
        S = diag(1./diag(sqrt(S)))*S*diag(1./diag(sqrt(S)));
        S(abs(S)<0.2) = 0;
        
        eigs = eig(S);
    end
    
    Sigmas(:,:,i) = S;
    Omegas(:,:,i) = inv(Sigmas(:,:,i));
end


%% Generating States and Task Data
%Generating Task Data. Column 1 is rest indicator, column 2 is task
%indicator
xDat = zeros(N,2,Subs);
xDat((N/2+1):N,2,:) = 1;
xDat(:,1,:) = 1-xDat(:,2,:);

% Generating transition matrix for each condition xDat = 0 & xDat = 1.
true_Q = zeros(num_states,num_states,2);

%transitions when xDat = 0
true_Q(:,:,1) = [0.98,0.02,0;...
    0.02,0.98,0;...
    0,0.5,0.5];

%transitions when xDat = 1
true_Q(:,:,2) = [0,0.5,0.5;...
    0,0.7,0.3;...
    0,0.02,0.98];

fprintf('Generating state sequence \n')
%Generating states
true_states = zeros(N,Subs);
for sub = 1:Subs
    true_states(:,sub) = states_from_exo(xDat(:,:,sub),num_states,true_Q);
end

%Reduce to single indicator variable
xDat = xDat(:,2,:);

%% Generating Data Frame
fprintf('Generating Data frame from multivariate gaussian \n')

Y = zeros(N,R,Subs);
mu = zeros(R,1);

for j = 1:Subs
    for i = 1:N
        Y(i,:,j) = mvnrnd(mu, Sigmas(:,:,true_states(i,j)));
    end
end

