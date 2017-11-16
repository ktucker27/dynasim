function [rho, L] = get_density(n, w, d, f, g, gamma)

L = full(master_matrix(n, w, d, f, g, gamma));
[~,~,v] = svd(L);

if rank(L) ~= size(L,1) - 1
    disp('WARNING: Liouvillian has multi-dimensional null space');
end

rho = wind(v(:,size(L,1)));
rho = rho/trace(rho);