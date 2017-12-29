function [tt, rho] = coh_two_time_corr(n, o, f, g, gamma, tau)

[L, ~, sp, sm] = coherent_matrix(n, o, f, g, gamma);
L = full(L);
[~,~,v] = svd(L);

if rank(L) ~= size(L,1) - 1
    disp('WARNING: Liouvillian has multi-dimensional null space');
end

rho = wind(v(:,size(L,1)));
rho = rho/trace(rho);

tt = two_time_calc(L, rho, tau, sp*2/n, sm*2/n);

end

function tt = two_time_calc(L, rho, tau, sp, sm)

tt = zeros(size(tau,1),1);

v = unwind(sm*rho);

dtau = tau(2) - tau(1);
if tau(1) ~= 0
    v = expm(tau(1)*L)*v;
end
del = expm(dtau*L);

for ti = 1:size(tau,1)
    m = sp*wind(v);
    
    tt(ti,1) = trace(m);
    
    %if mod(ti-1,5) == 0
    %    disp(['tau: ', num2str(tau(ti))]);
    %end
    
    v = del*v;
end

end