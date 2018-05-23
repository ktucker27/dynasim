function tt = exact_tt_grid_w(gvec, d, n, wvec, f, gamma, tau, savefile)

tt = cell(size(wvec,1), size(gvec,1));

for i=1:size(gvec)
    disp(['g = ', num2str(gvec(i))]);
    for j=1:size(wvec)
        disp(['w = ', num2str(wvec(j))]);
        
        % Since g is to be plotted on the x-axis of a color plot, it varies
        % with the columns of the return matrices
        tt{j,i} = [tau, two_time_corr(n, wvec(j,1), d, f, gvec(i,1), gamma, tau)];
        
        if nargin > 7
            save(savefile, 'gvec', 'dvec', 'n', 'w', 'f', 'gamma', 'tt');
        end
    end
end