function tt = two_time_corr(n, w, d, f, g, gamma, tau)

L = full(master_matrix(n, w, d, f, g, gamma));
[~,~,v] = svd(L);

if rank(L) ~= size(L,1) - 1
    disp('WARNING: Liouvillian has multi-dimensional null space');
end

rho = wind(v(:,2^(2*n)));
rho = rho/trace(rho);

tt = two_time_calc(L, rho, tau, n);

end

function tt = two_time_calc(L, rho, tau, n)

sp = [0,1;0,0];
sm = [0,0;1,0];

tt = zeros(size(tau,1),1);

for ti = 1:size(tau,1)
    sum = 0;
    for a = 1:n
        sap = get_sa(sp, n, a);
        for b=1:n
            sbm = get_sa(sm, n, b);
            sbmrhov = unwind(sbm*rho);
            v = expm(tau(ti)*L)*sbmrhov;
            m = sap*wind(v);
            sum = sum + trace(m);
        end
    end
    tt(ti,1) = sum/(n*n);
    
    if mod(ti,1000) == 0
        disp(['tau: ', num2str(tau(ti))]);
    end
end

end