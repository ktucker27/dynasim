function [csz, csp, csm] = get_full_coll_ops(n)

sz = [1,0;0,-1];
sp = [0,1;0,0];
sm = sp';

csz = zeros(2^n, 2^n);
csp = zeros(2^n, 2^n);
csm = zeros(2^n, 2^n);
for a=1:n
csz = csz + get_sa(sz, n, a);
csp = csp + get_sa(sp, n, a);
csm = csm + get_sa(sm, n, a);
end
csz = 0.5*csz;