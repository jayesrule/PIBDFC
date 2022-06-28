n_state = 2
n_covariate = 4

Xi_old = rand(n_state,n_state+n_covariate)
Z = [0,1;0,1;0,1;1,0;1,0;0,1;0,1;1,0;0,1;1,0]
Exo = randn([n_covariate,10])
a = 0
b = 10
S = 2
Y = [randn(4,1);5*randn(6,1)]
Omega = reshape([1,5],[1,1,2]);

Xi_new = Holsclaw_Param(Xi_old, Z, Exo, a, b)

[states, probs, t_prob] = Holsclaw_State(Xi_new, Z, Exo, Omega, Y)

