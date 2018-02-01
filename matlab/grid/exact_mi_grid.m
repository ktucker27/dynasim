function mi = exact_mi_grid(gvec, nvec, wvec, f, gamma, tau, rho0, savefile)

mi = cell(size(nvec,1), size(gvec,1));

for i=1:size(gvec)
    disp(['g = ', num2str(gvec(i))]);
    for j=1:size(nvec)
        disp(['n = ', num2str(nvec(j))]);
        
        n = nvec(j);
        L = full(master_matrix(n, wvec(j), zeros(n,1), f, gvec(i), gamma));
        
        % Since g is to be plotted on the x-axis of a color plot, it varies
        % with the columns of the return matrices
        mi{j,i} = [tau, mutual_info_vt(rho0,L,tau)];
        
        if nargin > 8
            save(savefile, 'gvec', 'nvec', 'wvec', 'f', 'gamma', 'rho0', 'mi');
        end
    end
end