function [map_states, sure] = mapstates(st_post)

N = size(st_post,2);
map_states = zeros(1,N);
sure = map_states;

for i = 1:N
    tbl = tabulate(st_post(:,i));
    sure(i) = max(tbl(:,3))/100;
    counts = tbl(:,2);
    [m,ind] = max(counts);
    map_states(i) = tbl(ind,1);
end

end
