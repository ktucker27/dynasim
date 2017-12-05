function tc = decay_time(tt, tol)

tcv = tt(abs(tt(:,2)) > tol*abs(tt(1,2)),1);
tc = tcv(end,1);

if tc > 0.95*tt(end,1)
    disp('WARNING: Function passed in to decay_time drops below tol within 5% of end time');
end