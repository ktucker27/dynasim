function [sz, sp, sm, szsz] = coherent_solve(n, o, f, g, gamma)

sz = zeros(size(o,1),1);
sp = zeros(size(o,1),1);
sm = zeros(size(o,1),1);
szsz = zeros(size(o,1),1);
for i=1:size(o,1)
    disp(['omega: ', num2str(o(i,1))]);
    [A, szm, spm, smm] = coherent_matrix_onsite(n, o(i), f, g, gamma);
    A = full(A);
    
    if rank(A) ~= size(A,1) - 1
        disp('WARNING: Liouvillian has multi-dimensional null space');
    end
    
    [~,~,v] = svd(A);
    rho = wind(v(:,size(A,1)));
    rho = rho/trace(rho);
    sz(i,1) = trace(rho*szm);
    sp(i,1) = trace(rho*spm);
    sm(i,1) = trace(rho*smm);
    szsz(i,1) = trace(rho*szm*szm);
end