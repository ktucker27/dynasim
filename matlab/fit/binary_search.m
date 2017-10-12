function idx = binary_search(v,x)

lidx = 1;
uidx = size(v,1);

if v(lidx) >= x
    idx = lidx;
    return;
end

if v(uidx) <= x
    idx = uidx;
    return;
end

while uidx - lidx > 1
    idx = floor((lidx + uidx)/2);
    if v(idx) < x
        lidx = idx;
    else
        uidx = idx;
    end
end

if abs(v(uidx) - x) < abs(v(lidx) - x)
    idx = uidx;
else
    idx = lidx;
end