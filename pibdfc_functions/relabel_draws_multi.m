%Function to relabel simulation results for easier comparison
function map = relabel_draws_multi(Omega_prev,st_prev, Omega_next,st_next, alpha)



%%
[p,~,k] =  size(Omega_prev);
[S,N] = size(st_prev);


parcor_prev = prec2parcor(Omega_prev);
parcor_next = prec2parcor(Omega_next);

%parcor = 1*(parcor~=0) ;


%%

%%
prev_laplace = parcor_prev;
next_laplace = parcor_next;
for i = 1:k
    prev_laplace(:,:,i) = prev_laplace(:,:,i)+diag(sum(prev_laplace(:,:,i)));
    next_laplace(:,:,i) = next_laplace(:,:,i)+diag(sum(next_laplace(:,:,i)));
end

mat_dist = zeros(k,k);
st_dist = mat_dist;
map = zeros(k,1);

for i = 1:k
    eig_est_i = eig(next_laplace(:,:,i));
    est_inds = find(st_next == i);
    
    l = length(est_inds);
    if l == 0
        l = 1;
    end
    
    for j = 1:k
        mat_dist(i,j) = sum((eig_est_i-eig(prev_laplace(:,:,j))).^2);
        
        
        true_ind_j = find(st_prev == j);

        st_dist(i,j) = 1 - sum(ismember(est_inds, true_ind_j))/l;
        
    end
end

mat_dist = mat_dist./max(max(mat_dist));

dist = alpha*mat_dist+(1-alpha)*st_dist;
%%
for st = 1:k
    [mins,a] = min(dist);
    [~,m] = min(mins);
   
    dist(:,m) = 10;
    dist(a(m),:) = 10;
    map(m) = a(m);
end

