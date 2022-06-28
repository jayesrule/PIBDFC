%Function to threshold the output of the Warnick model according to a
%nominal false detection rate.


function est_graph = bdfc_fdr_thresh(ryan,rate,plotting)

%% Arguements
%Inputs
%ryan: the output from the warnick model. Should be a structure
%rate: the nominal false detection rate with which to threshold
%plotting: logical indicator of whether to plot the thresholding value vs
%the fdr.

%Output
%est_graph: The thresholded adjacency matrix of the graph which holds the
%nominal fdr. 

%% Begin function

s = ryan.S;
ppi_edges = double.empty;

for i = 1:s
    ppi_edges = [ppi_edges; tril_vec(ryan.ppi_edges(:,:,i))];
end

failing = 1 - ppi_edges;
kappa = linspace(0,1,100);
fdr = zeros(size(kappa,1),size(kappa,2));

for i = 1:length(kappa)
    inds = ppi_edges >  kappa(i);
    fdr(i) = sum(failing(inds)/sum(inds));
end

if plotting
    
    plot(kappa,fdr)
    xlabel('threshold')
    ylabel('fdr')
    title('Choosing a threshold')
    
end

cands = kappa(fdr<rate);

threshold = cands(1);
est_graph = 1*(ryan.ppi_edges > threshold);

end
