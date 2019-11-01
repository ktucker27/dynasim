function [H, sz, sp, sm] = coherent_heff(n, o, f, g, gamma, gs)

if mod(n,2) ~= 0
    disp('ERROR - n must be even');
    return;
end

mvec = n/2:-1:-n/2;
sz = diag(mvec);
D = size(sz, 1);
sp = circshift(sqrt(diag((n/2 - mvec).*(n/2 + mvec + 1))),-1);
sm = circshift(sqrt(diag((n/2 + mvec).*(n/2 - mvec + 1))),1);
spm = sp*sm;
eyed = speye(D,D);

H = 0.5*o*(sp + sm) + 0.5*gamma*g*spm - 1i*0.5*(gamma*f*spm + gs*(n/2*eyed + sz));

