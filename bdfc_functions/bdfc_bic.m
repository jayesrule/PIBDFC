function bic = bdfc_bic(results)
p = size(results.Sig_save,1);
S = size(results.Sig_save,3);
T = size(results.ppi_HMM,2);
params = S*p*(p-1)/2 + 2*S^2 + T;

ll = 0;
Sig = mean(results.Sig_save,4);
for t=1:T
    max_state = find(results.ppi_HMM(:,t)==max(results.ppi_HMM(:,t)));
    ll = ll -2*logmvnpdf(results.data.TC(t,:),zeros(1,28),Sig(:,:,max_state));
end

bic = T*2*params+ll;


end

