function graph_est = fdr_est(kappa_post, fdr)
%%

%Input
%%kappa_post: a RxRxSxnsims array of the shrinkage coeficients
%%1/(1+lambda_ijs*tau_s)
%%fdr: real between 0 and 1, desired false discovery rate
%Output
%%graph_est: a RxRxS array of indicator variables for selection

kappa_thresh = squeeze(median(median(median(kappa_post,4),1),2));
%estimate of the inclusion probability
[R,~,S,~] = size(kappa_post); % R regions, S states


idx = tril(ones(R),-1) == 1; %indices to grab lower triangular of kappa

%Storage
fdr_final = zeros(S,1); %actual fdr for each state S
thresh_final = zeros(S,1); %threshold for ppi to select edge by state
graph_est = zeros(R,R,S);
ppi_est = zeros(R,R,S);
grid_size = 10000;

%Grid of thresholding variable
thresh_vec = linspace(0,1,grid_size);



%perform over states, different thresholds for each state?
for s = 1:S
    ppi =  mean(squeeze(kappa_post(:,:,s,:)<kappa_thresh(s)),3);  %estimate of shrinkage param for state s
    ppi_est(:,:,s) = ppi;
    fdr_star = zeros(1,grid_size); %storage to hold calculated fdr
    
    %iterate over grid of thresholding values
    for i = 1:grid_size
        d = 1*(ppi>= thresh_vec(i)); %decision for each element
        nu = 1-ppi; %probability the element is not selected
        d_vec = d(idx); %vectorize lower tri
        nu_vec = nu(idx); %vectorize lower tru
        
        %to avoid dividing by 0
        if sum(d_vec) == 0
            fdr_star(i) =  0;
        else
            %calculate fdr for the given threshold
            fdr_star(i) = sum(nu_vec.*d_vec)/sum(d_vec);
        end
    end
     
    id_keep = fdr_star ~= 0;
    fdr_star = fdr_star(id_keep);
    thresh_vec_s = thresh_vec(id_keep);
    %find the threshold the achieves the desired fdr
    [~,b] = min(abs(fdr-fdr_star)); %b: the index of the threshold that is closest to the desired fdr
    fdr_final(s) = fdr_star(b);
    thresh_final(s) = thresh_vec_s(b);
end


fprintf('Final FDR by state \n')
fdr_final

%return the selected graph in matrix form
for s = 1:S
    graph_est(:,:,s) = 1*(ppi_est(:,:,s) >= thresh_final(s));
end

%Put 1's along the diagonal for each state
graph_est(repmat(eye(R),[1,1,S])==1) = 1;

end