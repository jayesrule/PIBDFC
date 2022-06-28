%Function to evaluate multivariate normal liklihood
function llike = mvnloglike(X, mu, omega)
%% Input
%X nxp matrix of observations
%mu 1xp vector of means
%omega pxp precision (inverse covariance) matrix

%% Output
%llike scaler log likelihood

%% Begin Code

[n,p] = size(X);
like = zeros(n,1);
c = -p/2*log(2*pi)+.5*log(det(omega));
for i = 1:n
    d = X(i,:)'-mu';
    like(i) = c - 0.5*d'*omega*d;
end

llike = sum(like);
