function [es, ess, epp] = calc_ev(rho, sz, sp, sm)

sx = 0.5*(sp + sm);
sy = -1i*0.5*(sp - sm);

es = zeros(3, 1);
ess = zeros(3, 3);
es(1, 1) = trace(rho*sx);
es(2, 1) = trace(rho*sy);
es(3, 1) = trace(rho*sz);

epp = trace(rho*sp*sp);

ess(1, 1) = trace(rho*sx*sx);
ess(1, 2) = trace(rho*sx*sy);
ess(1, 3) = trace(rho*sx*sz);

ess(2, 1) = trace(rho*sy*sx);
ess(2, 2) = trace(rho*sy*sy);
ess(2, 3) = trace(rho*sy*sz);

ess(3, 1) = trace(rho*sz*sx);
ess(3, 2) = trace(rho*sz*sy);
ess(3, 3) = trace(rho*sz*sz);
