function tt = coll_two_time_corr_to(n, w, f, g, gamma, tau)

[cz, czz, cpm] = coll_steady_vals(n, w, f, gamma);

z0 = [cpm; 0.5*(1 + cz)];

A = coll_two_time_matrix_to(n, w, f, g, gamma);

dtau = tau(2) - tau(1);
dea = expm(dtau*A);

tt = zeros(size(tau,1),1);
z = z0;
for i=1:size(tau)
tt(i,1) = (1 - 1/n)*z(1,1) + (1/n)*z(2,1);
z = dea*z;
end