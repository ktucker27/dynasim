function tc = decay_time(tt, tol)

tcv = tt(abs(tt(:,2)) > tol*abs(tt(1,2)),1);
tc = tcv(end,1);