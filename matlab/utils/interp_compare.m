function D = interp_compare(x1, y1, x2, y2)

y1i = interp1(x1, y1, x2);        
D = max(abs(y2 - y1i));