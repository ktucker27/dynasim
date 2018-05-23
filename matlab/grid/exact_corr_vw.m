function G = exact_corr_vw(g, delta, n, wvec, f, gamma)

G = zeros(size(wvec,1), 1);

for i=1:size(wvec)
    disp(['w = ', num2str(wvec(i))]);
    [d,~] = even_lor(n, delta);
    L = full(master_matrix(n, wvec(i), d, f, g, gamma));
    [~,~,v] = svd(L);
    
    if rank(L) ~= size(L,1) - 1
        disp('WARNING: Liouvillian has multi-dimensional null space');
    end
    
    rho = wind(v(:,size(L,1)));
    rho = rho/trace(rho);
    
    G(i,1) = real(comp_corr(rho));
    %G(i,1) = real(comp_sz(rho));
end