function rho2 = partial_trace(rho, k)
% k indicates which particle to trace out (in [1,n]).  The bits of a row or
% column index are in particle major order.  We need to convert the 
% particle index to its corresponding bit (zero-based index used for bits)
n = log2(size(rho,1));
b = n - k;
rho2 = reduce(rho, b, 0) + reduce(rho, b, 1);
end

function B = reduce(A, k, val)
B = zeros(size(A,1)/2, size(A,2)/2);
for i=0:size(A,1)-1
    if get_bit(i,k) ~= val
        continue
    end
    
    i2 = remove_bit(i,k);
    for j=0:size(A,2)-1
        if get_bit(j,k) ~= val
            continue
        end
        
        j2 = remove_bit(j,k);
        
        B(i2+1, j2+1) = B(i2+1, j2+1) + A(i+1,j+1);
    end
end
end

function b = get_bit(idx, k)
b = mod(floor(idx*2^-k),2);
%disp(['get_bit(', num2str(idx), ',', num2str(k), ') = ', num2str(b)]);
end

function idx2 = remove_bit(idx,k)
upper = floor(idx*2^-(k+1))*2^k;
lower = idx - 2*upper - (2^k)*get_bit(idx,k);
idx2 = upper + lower;
%disp(['remove_bit(', num2str(idx), ',', num2str(k), ') = ', num2str(idx2)]);
end