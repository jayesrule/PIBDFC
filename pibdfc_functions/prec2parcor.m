%Function to compute the partial correlation matrix from the precision
%matrix
function parcor = prec2parcor(Omega)

[p,~,K] = size(Omega);

parcor = repmat(eye(p),[1,1,K]);

for k = 1:K
    for i = 2:p
        for j = 1:(i-1)
            parcor(i,j,k) = -Omega(i,j,k)/sqrt(Omega(i,i,k)*Omega(j,j,k));
            parcor(j,i,k) = parcor(i,j,k);
        end
    end
end

end
        
