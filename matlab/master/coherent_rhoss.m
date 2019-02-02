function rhoss = coherent_rhoss(o, f, g, smc, spc)

alpha = 1i*o/(f - 1i*g);

n = size(smc,1) - 1;

sp2 = alpha*eye(n+1,n+1) - spc;
sm2 = alpha'*eye(n+1,n+1) - smc;

rhoss = inv(sm2)*inv(sp2);

rhoss = rhoss/trace(rhoss);