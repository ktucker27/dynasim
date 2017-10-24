function v = unwind(rho)

n = size(rho,1);
v = zeros(n*n,1);

idx = 1;
for i=1:n
    for j = 1:n
        v(idx) = rho(i,j);
        idx = idx + 1;
    end
end