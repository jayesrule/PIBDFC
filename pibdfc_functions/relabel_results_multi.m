%Function to relabel simulation results for easier comparison
function [post_out] = relabel_results_multi(post_draws, Y, xDat, true_graph, true_states)



%%
kappa_post = post_draws.kappa_post;
parcor = fdr_est(kappa_post, 0.2);
[p,~,k] = size(parcor);


true_graph(repmat(eye(p),1,1,k)>0) = 0;
%parcor = 1*(parcor~=0) ;
alpha = 0.3;


%%
Omega_out = post_draws.Omega_post;
kappa_out = kappa_post;
rho_out = post_draws.rho_post;
eta_out = post_draws.eta_post;
Z_out = post_draws.Z_post;
sqrl_out = post_draws.sqrl_post;
st_post = post_draws.states_post;



%%

map_states = mapstates_mode(st_post);
true_laplace = true_graph;
est_laplace = parcor;
for i = 1:k
    true_laplace(:,:,i) = true_laplace(:,:,i)+diag(sum(true_laplace(:,:,i)));
    est_laplace(:,:,i) = est_laplace(:,:,i)+diag(sum(est_laplace(:,:,i)));
end

mat_dist = zeros(k,k);
st_dist = mat_dist;
map = zeros(k,1);

for i = 1:k
    eig_est_i = eig(est_laplace(:,:,i));
    
    for j = 1:k
        mat_dist(i,j) = sum((eig_est_i-eig(true_laplace(:,:,j))).^2);
        
        st_dist(i,j) = 1 - mean(map_states(true_states==j)==i);
        
    end
    
    
end

mat_dist = mat_dist./max(max(mat_dist));

dist = alpha*mat_dist+(1-alpha)*st_dist;

[~,a] = min(dist);
if length(unique(a))==k
    Omega_out = post_draws.Omega_post(:,:,a,:);
    kappa_out = kappa_post(:,:,a,:);
    rho_out = post_draws.rho_post(a,:,:,:);
    eta_out = post_draws.eta_post(a,:,:);
    for m = 1:k
        Z_out(m,:,:)= post_draws.Z_post(a(m),:,:);
        Z_out(:,m,:) = post_draws.Z_post(:,a(m),:);
        sqrl_out(m,:,:,:) = post_draws.sqrl_post(a(m),:,:,:);
        sqrl_out(:,m,:,:)= post_draws.sqrl_post(:,a(m),:,:);
    end
    
    map = a';
else
    
    %%
    for st = 1:k
        [mins,a] = min(dist);
        [~,m] = min(mins);
        
        Omega_out(:,:,m,:) = post_draws.Omega_post(:,:,a(m),:);
        kappa_out(:,:,m,:) = kappa_post(:,:,a(m),:);
        rho_out(m,:,:,:) = post_draws.rho_post(a(m),:,:,:);
        eta_out(m,:,:) = post_draws.eta_post(a(m),:,:);
        Z_out(m,:,:) = post_draws.Z_post(a(m),:,:);
        Z_out(:,m,:) = post_draws.Z_post(:,a(m),:);
        sqrl_out(m,:,:,:) = post_draws.sqrl_post(a(m),:,:,:);
        sqrl_out(:,m,:,:) = post_draws.sqrl_post(:,a(m),:,:);
        
        dist(:,m) = 10;
        dist(a(m),:) = 10;
        map(m) = a(m);
    end
    
end



inds = zeros([size(st_post),k]);
for i = 1:k
    inds(:,:,:,i) = st_post==i;
end

st_out = st_post;
for i = 1:k
    st_out(inds(:,:,:,i)==1) = find(map==i);
end

%lut = uint8([0; map; zeros(255 - length(map), 1)]);
%st_out = cast(intlut(uint8(st_post), lut), 'single') ;


post_out.Omega_post = Omega_out;
post_out.kappa_post = kappa_out;
post_out.rho_post = rho_out;
post_out.eta_post = eta_out;
post_out.sqrl_post = sqrl_out;
post_out.Z_post = Z_out;
post_out.states_post = st_out;

[[1:k]',map]
