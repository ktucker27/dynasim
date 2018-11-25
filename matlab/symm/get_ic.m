function c = get_ic(n)

d = (n+3)*(n+2)*(n+1)/6;

c = zeros(d,1);

col_idx = [0;0;0];
for i=1:d
     %if col_idx(2) == 0 && col_idx(3) == 0
         %c(i,1) = (1/2^n)*(-1)^col_idx(1);
         %c(i,1) = (1/2^n);
     %end
     if col_idx(1) == 0
         %c(i,1) = (1/2^n);
         c(i,1) = (1/2^n)*(-1)^(col_idx(2) + col_idx(3));
         %c(i,1) = (1/2^n)*(-1i)^col_idx(2)*(1i)^col_idx(3);
         %c(i,1) = (1/2^n)*(1i)^col_idx(2)*(-1i)^col_idx(3);
     end
     col_idx = inc_idx(col_idx, n);
end