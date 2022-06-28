%function to generate a block diagonal matrix of size p

function [Sigma] = Generateblkdiag(ncliques,cliquesize,theta)

p = ncliques * cliquesize;
sizes =  repelem(cliquesize, ncliques);
Sigma = zeros(p);
lims = [0, cumsum(sizes)];
for i = 1:ncliques 
        ii  = (lims(i) + 1):lims(i + 1);
        signs = 2 * (normrnd(0,1,sizes(i) * (sizes(i) - 1)/2,1) < 0) - 1;
        d = Sigma(ii, ii); 
        ut = triu(ones(size(d)),1)== 1;
        d(ut) =  signs * theta;
        Sigma(ii, ii) = d;
end

Sigma = Sigma + Sigma';
eigen = eig(Sigma);
shift = (max(eigen) - p * min(eigen))/(p - 1);
    %cat("Shifting eigenvalues by ", shift, fill = T)
Sigma(eye(p)==1) = diag(Sigma) + shift;
    %A <- eig$vect %*% diag(sqrt(eig$val + shift)) %*% t(eig$vect)
    %list(Sigma = Sigma, A = A, shift = shift)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    