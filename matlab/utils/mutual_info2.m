function info = mutual_info2(pms, zs, zzs, ps, zps, pps)

info = zeros(size(pms,1),1);
for j=1:size(pms,1)
    %info(j) = single_info(zs(j), 0) + ...
    %          single_info(zs(j), 0) - ...
    %          double_info(0, 0, zs(j), zs(j), 0, 0, 0, pms(j), zzs(j));
          
    info(j) = single_info(zs(j), ps(j)) + ...
                    single_info(zs(j), ps(j)) - ...
                    double_info(ps(j), ps(j), zs(j), zs(j), zps(j), zps(j), pps(j), pms(j), zzs(j));
end
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
if(abs(sum(l)-1.0) > 1.0e-10 || max(max(abs(A - A'))) > 1.0e-10 || min(l >= -1.0e-4) < 1)
sum(l)
A - A'
l >= 0
l
end
end