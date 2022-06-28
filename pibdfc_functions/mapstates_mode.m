function [map_states] = mapstates_mode(st_post)

NT = size(st_post,1);
nsubs = size(st_post,2);
map_states = zeros(NT, nsubs);

for sub = 1:nsubs
    for i = 1:NT
        tbl = tabulate(squeeze(st_post(i,sub,:)));
        counts = tbl(:,2);
        [~,ind] = max(counts);
        map_states(i,sub) = tbl(ind,1);
    end
end

end
