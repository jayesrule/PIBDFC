function out = tril_vec(A)
%Function to out put the lower triangular of a matrix as a vector
%// Mask of lower triangular elements
mask = tril(true(size(A)),-1);

%// Use mask to select lower triangular elements from input array
out = A(mask);
end
