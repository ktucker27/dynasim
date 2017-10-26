function v = unwind(rho)

n = size(rho,1);
v = reshape(transpose(rho), n*n, 1);