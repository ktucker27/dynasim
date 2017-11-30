function tt = coll_two_time_corr(n, z0, A, tau)

dtau = tau(2) - tau(1);
dea = expm(dtau*A);

tt = zeros(size(tau,1),1);
z = z0;
for i=1:size(tau)
tt(i,1) = (1 - 1/n)*z(1,1) + (1/n)*z(2,1);
z = dea*z;
end