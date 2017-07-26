function [xp,yp] = peak_picker(x,y)

di = 25;
xp = [];
yp = [];

while max(abs(y)) > 100*mean(abs(y))
    [~,i] = max(abs(y));
    xp = [xp;x(i)];
    yp = [yp;y(i)];
    y(i-di:i+di) = 0;
end