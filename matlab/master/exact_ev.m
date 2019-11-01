function ev = exact_ev(rho, L, op, t)

dt = 0;
if size(t,1) > 1
    dt = t(2) - t(1);
end

del = expm(dt*L);

if t(1) > 0
    del0 = expm(t(1)*L);
    rho = wind(del0*unwind(rho));
end

ev = zeros(size(t,1),1);
for i=1:size(t,1)
    ev(i,1) = trace(rho*op);
    rho = wind(del*unwind(rho));
end