%function to half vectorize a matrix: A
function vec_h = half_vec(A)

vec_h = A(tril(true(size(A))));

end