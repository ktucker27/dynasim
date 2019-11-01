function [es, ess, epp, psiss, norms] = get_all_ev_vec(psi, H, sz, sp, sm, t0, dt, t)

sx = 0.5*(sp + sm);
sy = -1i*0.5*(sp - sm);

del = expm(-1i*dt*H);

if t0 ~= 0
    psi = expm(-1i*t0*H)*psi;
end

numt = round((t - t0)/dt + 1);

%psin = psi;

epp = zeros(numt,1);
es = zeros(3,numt);
ess = zeros(3, 3, numt);
norms = zeros(numt,1);
for i=1:numt
    psin = psi/norm(psi,2);
    
    es(1, i) = psin'*sx*psin;
    es(2, i) = psin'*sy*psin;
    es(3, i) = psin'*sz*psin;
    
    epp(i,1) = psin'*sp*sp*psin;
    
    ess(1, 1, i) = psin'*sx*sx*psin;
    ess(1, 2, i) = psin'*sx*sy*psin;
    ess(1, 3, i) = psin'*sx*sz*psin;
    
    ess(2, 1, i) = psin'*sy*sx*psin;
    ess(2, 2, i) = psin'*sy*sy*psin;
    ess(2, 3, i) = psin'*sy*sz*psin;
    
    ess(3, 1, i) = psin'*sz*sx*psin;
    ess(3, 2, i) = psin'*sz*sy*psin;
    ess(3, 3, i) = psin'*sz*sz*psin;
    
    psi = del*psi;
    
    norms(i,1) = norm(psi,2);
    
    if norm(psi,2) < .2
        psi = psi/norm(psi,2);
    end
    
    tval = t0 + i*dt;
    %psin = expm(-1i*tval*H)*psi;
    %psin = psin/norm(psin,2);
    %disp(tval);
end

psiss = psi;