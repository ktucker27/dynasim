function A = master_matrix(n, w, d, f, g, gamma)

sz = [1,0;0,-1];
sp = [0,1;0,0];
sm = sp';
smp = sm*sp;
D = 2^n;
eyed = speye(D,D);
A = sparse(D^2, D^2);
for a = 1:n
    sap = get_sa(sp,n,a);
    A = A + 1i*0.5*d(a)*(kron(get_sa(sz,n,a), eyed) - kron(eyed, get_sa(sz,n,a)));
    A = A - 0.5*w*(kron(get_sa(smp,n,a), eyed) + kron(eyed, get_sa(smp,n,a)) - 2*kron(sap, sap));
    for b = 1:n
        sapsbm = get_sa(sp,n,a)*get_sa(sm,n,b);
        
        A = A - 0.5*gamma*f*(kron(sapsbm, eyed) + kron(eyed, sapsbm') - 2*(kron(get_sa(sm,n,b), sap')));
        
        if a ~= b
            A = A - 1i*0.5*gamma*g*(kron(sapsbm, eyed) - kron(eyed, sapsbm'));
        end
    end
end