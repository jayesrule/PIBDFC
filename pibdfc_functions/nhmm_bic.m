function bic = nhmm_bic(Y, Omega_est,map_states)

n = size(Y,1);
[p,~,k] = size(Omega_est);
mu = zeros(1,p);
nparams = k*p*(p+1)/2+k*(k+1);

for i = 1:n
    llike = mvnloglike(Y(i,:),mu,Omega_est(:,:,map_states(i)));
end

llike = sum(sum(llike));

bic = log(n)*nparams-2*llike;



