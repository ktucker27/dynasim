function [es, ess, epp, rhoss] = coherent_ev(rho, L, sz, sp, sm, t0, dt, t)

sx = 0.5*(sp + sm);
sy = -1i*0.5*(sp - sm);

del = expm(dt*L);

if t0 ~= 0
    rho = wind(expm(t0*L)*unwind(rho));
end

numt = (t - t0)/dt + 1;

epp = zeros(numt,1);
es = zeros(3,numt);
ess = zeros(3, 3, numt);
for i=1:numt
    es(1, i) = trace(rho*sx);
    es(2, i) = trace(rho*sy);
    es(3, i) = trace(rho*sz);
    
    epp(i,1) = trace(rho*sp*sp);
    
    ess(1, 1, i) = trace(rho*sx*sx);
    ess(1, 2, i) = trace(rho*sx*sy);
    ess(1, 3, i) = trace(rho*sx*sz);
    
    ess(2, 1, i) = trace(rho*sy*sx);
    ess(2, 2, i) = trace(rho*sy*sy);
    ess(2, 3, i) = trace(rho*sy*sz);
    
    ess(3, 1, i) = trace(rho*sz*sx);
    ess(3, 2, i) = trace(rho*sz*sy);
    ess(3, 3, i) = trace(rho*sz*sz);
    
    rho = wind(del*unwind(rho));
end

rhoss = rho;