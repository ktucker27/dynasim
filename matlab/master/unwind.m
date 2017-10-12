function rho = unwind(v)

n = sqrt(size(v,1));
rho = zeros(n,n);

idx = 1;
for i=1:n
    for j = 1:n
        rho(i,j) = v(idx);
        idx = idx + 1;
    end
end