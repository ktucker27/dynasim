function [b,ok,diff] = exp_fit(t,z,tol,ts)

if nargin < 4
    ts = 1;
end

tsi = binary_search(t,ts);
t = t(tsi:end);
z = z(tsi:end);

[x,y] = plot_fft(t,z,0);
[xp,yp] = peak_picker(x,y);

b0 = zeros(2*size(xp,1),1);
idx = 1;
for i=1:size(xp,1)
    b0(idx) = yp(i)/size(z,1);
    b0(idx+1) = xp(i);
    idx = idx + 2;
end

fun = @(b,t)exp_comb(b,t);

try
    ops.MaxIter = 1000;
    ops.TolFun = 1e-15;
    b = nlinfit(t, z, fun, b0, ops);
catch ex
    fprintf('ERROR - Nonlinear fit failed:\n%s\n', getReport(ex));
    ok = 0;
    diff = Inf;
    return;
end
diff = verify_fit(t,z,b,fun);

if(diff >= tol)
    fprintf('WARNING - Fit failed with diff %g\n', diff);
    ok = 0;
else
    ok = 1;
end
end

function diff = verify_fit(t,z,b,fun)
y = fun(b,t);
diff = max(abs(y - z));
end