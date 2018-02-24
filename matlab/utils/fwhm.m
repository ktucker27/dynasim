function B = fwhm(x, y, xprec)

if x(2) - x(1) > xprec
    x2 = (x(1):xprec:x(size(x,1)))';
    y = interp1(x, y, x2);
    x = x2;
end

ind = (abs(y) > 0.5*max(abs(y)));
idx1 = find(ind,1);
idx2 = find(ind,1,'last');
B = x(idx2) - x(idx1);
