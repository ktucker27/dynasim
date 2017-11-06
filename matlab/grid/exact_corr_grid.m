function G = exact_corr_grid(gvec, dvec, n, w, f, gamma)

G = zeros(size(gvec,1), size(dvec,1));

for i=1:size(gvec)
    for j=1:size(dvec)
        [d,~] = even_lor(n, dvec(j));
        L = full(master_matrix(n, w, d, f, gvec(i), gamma));
        [~,~,v] = svd(L);
        
        if rank(L) ~= size(L,1) - 1
            disp('WARNING: Liouvillian has multi-dimensional null space');
        end
        
        rho = wind(v(:,size(L,1)));
        rho = rho/trace(rho);
        
        G(j,i) = real(comp_corr(rho));
    end
end