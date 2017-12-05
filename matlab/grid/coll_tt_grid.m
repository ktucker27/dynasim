function tt = coll_tt_grid(gvec, n, w, f, gamma, tau, to)

tt = cell(1, size(gvec,1));

for i=1:size(gvec)
    disp(['g = ', num2str(gvec(i))]);
    for j=1:size(w,1)
        disp(['w = ', num2str(w(j))]);
        if nargin > 6 && to ~= 0
            tt{j,i} = [tau, coll_two_time_corr_to(n, w(j), f, gvec(i), gamma, tau)];
        else
            tt{j,i} = [tau, coll_two_time_corr(n, w(j), f, gvec(i), gamma, tau)];
        end
    end
end