function [es, ess] = mcwf_ev(n, rho)

numt = size(rho,3);

es = zeros(3, numt);
ess = zeros(3, 3, numt);

dim = size(rho,1);
sz = zeros(dim, dim);

idx = 1;
for J = n/2:-1:-n/2
    for M = J:-1:-J
        sz(idx, idx) = M;
        idx = idx + 1;
    end
end

for ti = 1:numt
    es(3,ti) = trace(sz*rho(:,:,ti));
end