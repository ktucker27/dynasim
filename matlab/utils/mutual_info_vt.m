function info = mutual_info_vt(rho, L, t)

n = log2(size(rho,1));

info = zeros(size(t,1),1);

dt = t(2) - t(1);

del = expm(dt*L);
for i=1:size(t,1)
    % Trace down to two particle density operator
    rhoAB = rho;
    for j=n:-1:3
        rhoAB = partial_trace(rhoAB, j);
    end
    l = eig(rhoAB);
    sab = -sum(l.*log(l));
    
    % Trace down to single particle density operator
    rhoA = partial_trace(rhoAB,2);
    l = eig(rhoA);
    sa = -sum(l.*log(l));
    
    rhoB = partial_trace(rhoAB,1);
    l = eig(rhoB);
    sb = -sum(l.*log(l));
    
    info(i,1) = sa + sb - sab;
    
    rho = wind(del*unwind(rho));
end
end
