function tc = decay_time(tt, tol)

tcv = tt(abs(tt(:,2)) > tol,1);
tc = tcv(end,1);