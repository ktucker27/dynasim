function tt = two_time_corr(n, w, d, f, g, gamma, tau)

L = full(master_matrix(n, w, d, f, g, gamma));
[~,~,v] = svd(L);

if rank(L) ~= size(L,1) - 1
    disp('WARNING: Liouvillian has multi-dimensional null space');
end

rho = wind(v(:,size(L,1)));
rho = rho/trace(rho);

tt = two_time_calc(L, rho, tau, n);

end

function tt = two_time_calc(L, rho, tau, n)

sp = [0,1;0,0];
sm = [0,0;1,0];

tt = zeros(size(tau,1),1);

sum1 = zeros(2^n,2^n);
sum2 = zeros(2^n,2^n);
for a = 1:n
    sum1 = sum1 + get_sa(sp, n, a);
    sum2 = sum2 + get_sa(sm, n, a);
end

v = unwind(sum2*rho);

dtau = tau(2) - tau(1);
if tau(1) ~= 0
    v = expm(tau(1)*L)*v;
end
del = expm(dtau*L);

for ti = 1:size(tau,1)
    m = sum1*wind(v);
    
    tt(ti,1) = trace(m)/(n*n);
    
    %if mod(ti-1,5) == 0
    %    disp(['tau: ', num2str(tau(ti))]);
    %end
    
    v = del*v;
end

end