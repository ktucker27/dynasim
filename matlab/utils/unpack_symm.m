function [tvec, es, ess, es2, ess2] = unpack_symm(B)

tvec = B(:,1);
nt = size(tvec,1);

sq_disp = 21;

es = [B(:,2)';B(:,3)';B(:,4)'];
if size(B,2) > 22
    es2 = [B(:,2 + sq_disp)';B(:,3 + sq_disp)';B(:,4 + sq_disp)'];
else
    es2 = zeros(size(es));
end

ess = zeros(3,3,nt);
ess2 = zeros(3,3,nt);
for i=1:nt
    idx = 5;
    for j=1:3
        for k=1:3
            ess(j,k,i) = B(i,idx) + 1i*B(i,idx+1);
            if size(B,2) > 22
                ess2(j,k,i) = B(i,idx + sq_disp) + 1i*B(i,idx+1 + sq_disp);
            end
            idx = idx + 2;
        end
    end
end