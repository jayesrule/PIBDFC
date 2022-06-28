function [ mrf_C ] = calc_mrf_C(Theta, nu)
% Calculate normalizing constant for MRF prior for specific edge (i, j)

mrf_C = 0;
K = size(Theta, 1);

% Generate matrix with all possible binary vectors as the rows
% From http://www.mathworks.com/matlabcentral/newsreader/view_thread/121631
D = (0:2^K-1)';
B = rem(floor(D * pow2(-(K - 1):0)), 2);

for i = 1:size(B, 1)
    g_ij = B(i, :);
    mrf_C = mrf_C + exp(nu * sum(g_ij) + (g_ij) * Theta * (g_ij)');
end


end

