function [sq2mat, omat] = sq2_ss_grid(nvec, oratio, f, g)

sq2mat = zeros(size(nvec,1), size(oratio,1));
omat = zeros(size(nvec,1), size(oratio,1));
for i=1:size(nvec)
    n = nvec(i);
    oc = (n/2)*sqrt(f^2 + 4*(g/2)^2);
    for j=1:size(oratio)
        o = oratio(j)*oc;
        disp(['Processing: n = ', num2str(n), ', o/oc = ', num2str(oratio(j))]);
        sq2mat(i,j) = get_ss_sq2(n,o,f,g);
        omat(i,j) = o;
    end
end