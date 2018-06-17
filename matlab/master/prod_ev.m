function ev = prod_ev(n1, n2, ess)

ev = 0;
for i=1:3
    for j=1:3
        ev = ev + n1(i)*n2(j)*ess(i,j);
    end
end
