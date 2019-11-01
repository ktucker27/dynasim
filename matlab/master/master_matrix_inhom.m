function [A, szc, spc, smc] = master_matrix_inhom(n, w, o, d, fmat, gmat, gvec, gamma, gel)

sz = [1,0;0,-1];
sp = [0,1;0,0];
sm = sp';
smp = sm*sp;
D = 2^n;
eyed = speye(D,D);
A = sparse(D^2, D^2);
szc = zeros(D,D);
spc = zeros(D,D);
smc = zeros(D,D);
for a = 1:n
    sap = get_sa(sp,n,a);
    sam = get_sa(sm,n,a);
    saz = get_sa(sz,n,a);
    
    szc = szc + (1/2)*saz;
    spc = spc + sap;
    smc = smc + sam;
    
    % Disorder
    A = A + 1i*0.5*d(a)*(kron(saz, eyed) - kron(eyed, saz));
    
    % Incoherent pumping
    A = A - 0.5*w*(kron(get_sa(smp,n,a), eyed) + kron(eyed, get_sa(smp,n,a)) - 2*kron(sap, sap));
    
    for b = 1:n
        sapsbm = get_sa(sp,n,a)*get_sa(sm,n,b);
        
        f = fmat(a,b);
        g = gmat(a,b);
        
        % Inelastic interactions
        A = A - 0.5*gamma*f*(kron(sapsbm, eyed) + kron(eyed, sapsbm') - 2*(kron(get_sa(sm,n,b), sap')));
        
        % Elastic interactions
        A = A - 1i*0.5*gamma*g*(kron(sapsbm, eyed) - kron(eyed, sapsbm'));
    end
    
    % Coherent pumping
    A = A - 1i*0.5*gvec(a)*o*(kron(sap, eyed) + kron(sam, eyed) - kron(eyed, sam) - kron(eyed, sap));
    
    % Spontaneous emission
    A = A - 0.5*gel*(kron(eyed, eyed) - kron(saz, saz));
end