function [tvec, es, ess] = unpack_symm(B)

tvec = B(:,1);
nt = size(tvec,1);

es = [B(:,2)';B(:,3)';B(:,4)'];

ess = zeros(3,3,nt);
for i=1:nt
    idx = 5;
    for j=1:3
        for k=1:3
            ess(j,k,i) = B(i,idx) + 1i*B(i,idx+1);
            idx = idx + 2;
        end
    end
end