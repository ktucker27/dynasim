function G = exact_corr_vw(g, delta, n, wvec, f, gamma)

G = zeros(size(wvec,1), 1);

for i=1:size(wvec)
    disp(['w = ', num2str(wvec(i))]);
    [d,~] = even_lor(n, delta);
    L = full(master_matrix(n, wvec(i), 0, d, f, f, 0, g, gamma, 0));
    [~,~,v] = svd(L);
    
    if rank(L) ~= size(L,1) - 1
        disp('WARNING: Liouvillian has multi-dimensional null space');
    end
    
    rho = wind(v(:,size(L,1)));
    rho = rho/trace(rho);
    
    spsm = real(comp_corr(rho));
    sz = real(comp_sz(rho));
    
    G(i,1) = (1/n^2)*(n*(n-1)*spsm + (n/2)*(1 + sz));
end