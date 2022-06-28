%Function to relabel simulation results for easier comparison
function [st_out, map] = relabel_results2(map_states, true_states)



%%
k = max(max(true_states));
N = max(size(map_states));


frob_dist = zeros(k,k);
map = zeros(k,1);

for i = 1:k
    for j = 1:k
        matched = (map_states==i) & (true_states==j);
        frob_dist(i,j) = sum(sum(matched==0));
    end
end

dist = frob_dist/N;
%%
% for st = 1:k
%     [mins,a] = min(dist);
%     [~,m] = min(mins);
%     dist(:,m) = 10;
%     dist(a(m),:) = 10;
%     map(a(m)) = m;
% end

[~,a] = min(dist);
map = a';
inds = zeros(N,k);
for i = 1:k
    inds(:,i) = map_states==i;
end

st_out = map_states;
for i = 1:k
    st_out(inds(:,i)==1) = find(map==i);
end

[[1:k]',map]
