function ev = exact_ev(rho, L, op, t)

dt = t(2) - t(1);

del = expm(dt*L);

ev = zeros(size(t,1),1);
for i=1:size(t,1)
    ev(i,1) = trace(rho*op);
    rho = wind(del*unwind(rho));
end