function rho = wind(v)

n = sqrt(size(v,1));
rho = transpose(reshape(v,n,n));
