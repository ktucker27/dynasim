function [es, ess, epp] = compact_ev(n, L, c, t0, dt, t)

del = expm(dt*L);

if t0 ~= 0
    c = expm(t0*L)*c;
end

idx100 = get_idx([1;0;0], n);
idx010 = get_idx([0;1;0], n);
idx001 = get_idx([0;0;1], n);
idx200 = get_idx([2;0;0], n);
idx020 = get_idx([0;2;0], n);
idx002 = get_idx([0;0;2], n);
idx110 = get_idx([1;1;0], n);
idx101 = get_idx([1;0;1], n);
idx011 = get_idx([0;1;1], n);

nc2 = nchoosek(n,2);

numt = (t - t0)/dt + 1;
epp = zeros(numt, 1);
es = zeros(3, numt);
ess = zeros(3, 3, numt);
for i=1:numt
    es(1,i) = (n/2)*2^(n-1)*(c(idx001) + c(idx010));
    es(2,i) = 1i*(n/2)*2^(n-1)*(-c(idx001) + c(idx010));
    es(3,i) = (n/2)*2^n*c(idx100);
    
    epp(i,1) = 2^(n-1)*nc2*c(idx002);
    
    ess(1, 1, i) = 0.25*(2^(n-1)*nc2*(c(idx020) + c(idx002)) ...
                 + n*(1 + (n-1)*2^(n-1)*c(idx011)));
             
    ess(1, 2, i) = 1i*2^(n-3)*nc2*(c(idx020) - c(idx002)) + 1i*n*2^(n-2)*c(idx100);
    
    ess(2, 1, i) = 1i*2^(n-3)*nc2*(c(idx020) - c(idx002)) - 1i*n*2^(n-2)*c(idx100);
             
    ess(1, 3, i) = 2^(n-3)*n*((n-1)*(c(idx110) + c(idx101)) - c(idx001) + c(idx010));
    
    ess(3, 1, i) = 2^(n-3)*n*((n-1)*(c(idx110) + c(idx101)) + c(idx001) - c(idx010));
             
    ess(2, 2, i) = -0.25*(2^(n-1)*nc2*(c(idx020) + c(idx002)) ...
                 - n*(1 + (n-1)*2^(n-1)*c(idx011)));
             
    ess(2, 3, i) = -1i*2^(n-3)*n*((n-1)*(c(idx101) - c(idx110)) - c(idx001) - c(idx010));
    
    ess(3, 2, i) = -1i*2^(n-3)*n*((n-1)*(c(idx101) - c(idx110)) + c(idx001) + c(idx010));
             
    ess(3, 3, i) = n/4 + 2^(n-1)*nc2*c(idx200);
    
    c = del*c;
end
end

function idx = get_idx(idx_vec, n)
d = (n+3)*(n+2)*(n+1)/6;
comp_vec = [0;0;0];
idx = 1;
for i=1:d
    if min(idx_vec == comp_vec) == 1
        idx = i;
        break;
    end
    
    comp_vec = inc_idx(comp_vec, n);
end
end