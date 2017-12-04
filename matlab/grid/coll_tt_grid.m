function tt = coll_tt_grid(gvec, n, w, f, gamma, tau, to)

tt = cell(1, size(gvec,1));

for i=1:size(gvec)
    disp(['g = ', num2str(gvec(i))]);
    
    if nargin > 6 && to ~= 0
        tt{1,i} = [tau, coll_two_time_corr_to(n, w, f, gvec(i), gamma, tau)];
    else
        tt{1,i} = [tau, coll_two_time_corr(n, w, f, gvec(i), gamma, tau)];
    end
end