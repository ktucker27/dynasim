function idx = inc_idx(idx, n)
if min(idx == -1) == 1
    return;
end

idx(1) = idx(1) + 1;
if(sum(idx) > n)
    if(idx(1) ~= 1)
        idx(1) = 0;
        idx(2) = idx(2) + 1;
    else
        idx(1) = 0;
        idx(2) = 0;
        idx(3) = idx(3) + 1;
        if(idx(3) > n)
            idx = [-1;-1;-1];
        end
    end
end
end