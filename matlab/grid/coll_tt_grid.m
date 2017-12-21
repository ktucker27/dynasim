function tt = coll_tt_grid(gvec, n, w, f, gamma, tau, to)

tt = cell(1, size(gvec,1));

if size(w,1) > 1 && size(n,1) == 1
    n = n*ones(size(w,1),1);
end

for i=1:size(gvec)
    disp(['g = ', num2str(gvec(i))]);
    for j=1:size(w,1)
        disp(['w = ', num2str(w(j))]);
        if nargin > 6 && to ~= 0
            tt{j,i} = [tau, coll_two_time_corr_to(n(j), w(j), f, gvec(i), gamma, tau)];
        else
            tt{j,i} = [tau, coll_two_time_corr(n(j), w(j), f, gvec(i), gamma, tau)];
        end
    end
end