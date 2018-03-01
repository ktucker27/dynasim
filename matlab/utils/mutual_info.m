function [avg_info, info] = mutual_info(filepath, n)

B = dlmread(filepath);
z = B(1:2:end,1) + 1i*B(2:2:end,1);
[~, ps, zps, zs, pms, zzs, pps] = z_to_vars(z, n);

info = zeros(n,n);
for i=1:n
    for j=1:n
        if(i == j)
            continue;
        end
        info(i,j) = single_info(zs(i), ps(i)) + ...
                    single_info(zs(j), ps(j)) - ...
                    double_info(ps(i), ps(j), zs(i), zs(j), zps(i,j), zps(j,i), pps(i,j), pms(i,j), zzs(i,j));
    end
end

avg_info = sum(sum(info))/(n*n - n);
end

function info = single_info0(sz, sp)
A = [0.5 - sz, sp;sp', 0.5 + sz];
l = eig(A);
info = -1.0*sum(l.*log(l));
end

function info = double_info0(a_p, b_p, a_z, b_z, ab_zp, ba_zp, pp, pm, zz)
A = [0.0, 0.5*(b_p - 2*ab_zp), 0.5*(a_p - 2*ba_zp), pp; ...
     0.0, 0.0, pm, 0.5*a_p + ba_zp; ...
     0.0, 0.0, 0.0, 0.5*b_p + ab_zp; ...
     0.0, 0.0, 0.0, 0.0];
A = A + A';
A(1,1) = 0.25*(1 - 2*a_z - 2*b_z) + zz;
A(2,2) = 0.25*(1 - 2*a_z + 2*b_z) - zz;
A(3,3) = 0.25*(1 + 2*a_z - 2*b_z) - zz;
A(4,4) = 0.25*(1 + 2*a_z + 2*b_z) + zz;
l = eig(A);
info = -1.0*sum(l.*log(l));
end

function info = single_info(sz, sp)
A = [0.5*(1 + sz), sp';sp, 0.5*(1 - sz)];
l = eig(A);
info = -1.0*sum(l.*log(l));
end

function info = double_info(a_p, b_p, a_z, b_z, ab_zp, ba_zp, pp, pm, zz)
A = [0.0, 2*b_p' + 2*ab_zp', 2*a_p' + 2*ba_zp', 4*pp'; ...
     0.0, 0.0, 4*pm', 2*a_p' - 2*ba_zp'; ...
     0.0, 0.0, 0.0, 2*b_p' - 2*ab_zp'; ...
     0.0, 0.0, 0.0, 0.0];
A = A + A';
A(1,1) = 1 + a_z + b_z + zz;
A(2,2) = 1 + a_z - b_z - zz;
A(3,3) = 1 - a_z + b_z - zz;
A(4,4) = 1 - a_z - b_z + zz;
A = 0.25*A;
l = eig(A);
info = -1.0*sum(l.*log(l));
if(abs(sum(l)-1.0) > 1.0e-10 || max(max(abs(A - A'))) > 1.0e-10 || min(l >= 0) < 1)
sum(l)
A - A'
l >= 0
end
end