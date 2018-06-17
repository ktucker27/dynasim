function [A, sz, sp, sm] = coherent_matrix_onsite(n, o, f, g, gamma)

mvec = n/2:-1:-n/2;
sz = diag(mvec);
D = size(sz, 1);
sp = circshift(sqrt(diag((n/2 - mvec).*(n/2 + mvec + 1))),-1);
sm = circshift(sqrt(diag((n/2 + mvec).*(n/2 - mvec + 1))),1);
spm = sp*sm;
eyed = speye(D,D);
A = sparse(D^2, D^2);

A = A - 1i*0.5*o*(kron(sp, eyed) + kron(sm, eyed) - kron(eyed, sm) - kron(eyed, sp));
A = A - 1i*0.5*gamma*g*(kron(spm, eyed) - kron(eyed, spm));
A = A - 0.5*gamma*f*(kron(spm, eyed) + kron(eyed, spm) - 2*kron(sm, eyed)*kron(eyed, sm));
