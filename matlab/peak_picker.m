function [xp,yp] = peak_picker(x,y)

di = 20;
xp = [];
yp = [];

[v,i] = max(abs(y));
start_max = v;

%while max(abs(y)) > 100*mean(abs(y))
while v > 0.03*start_max
    xp = [xp;x(i)];
    yp = [yp;y(i)];
    y(i-di:i+di) = 0;
    [v,i] = max(abs(y));
end