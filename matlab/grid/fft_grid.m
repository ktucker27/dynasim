function [V, M, B, T, G, Dx, Dy] = fft_grid(tt, tol)

if nargin == 1
    tol = 0.001;
end

V = zeros(size(tt,1), size(tt,2));
M = zeros(size(tt,1), size(tt,2));
B = zeros(size(tt,1), size(tt,2));
T = zeros(size(tt,1), size(tt,2));
Dx = zeros(size(tt,1), size(tt,2));
Dy = zeros(size(tt,1), size(tt,2));

for i=1:size(tt,1)
    for j=1:size(tt,2)
        %[x,y] = plot_fft(tt{i,j}(:,1), tt{i,j}(:,2), 0);
        z1 = tt{i,j}(:,2);
        t1 = tt{i,j}(:,1);
        dt = t1(2,1) - t1(1,1);
        t = (t1(1,1):dt:t1(end,1))';
        z = interp1(t1, z1, t);
        %t = t1;
        %z = z1;
        %z = [z;zeros(3*size(z,1),1)];
        %t = [t;zeros(3*size(t,1),1)];
        [x,y] = plot_fft(t, z, 0);
        [x2,y2] = plot_fft_shift(t, z, 0);
        Dx(i,j) = max(abs(x - x2));
        Dy(i,j) = max(abs(y - y2));
        
        y = (abs(y) > 0.07*max(abs(y))).*y;
        
        [V(i,j), M(i,j), B(i,j)] = disc_dist(x,y);
        %[V(i,j), M(i,j), B(i,j)] = smooth_max(x,y);
        B(i,j) = fwhm(x, y, 0.001);
        T(i,j) = decay_time(tt{i,j}, tol);
        %M(i,j) = tt{i,j}(1,2);
    end
end

G = -log(tol)./T;

end

function [v,m,b] = smooth_max(x,y)
    %ay = smooth(abs(y), 3);
    ay = abs(y);
        
    [maxval, maxidx] = max(ay);
    
    v = x(maxidx);
    m = maxval;
    b = whm(x,ay);
end

function [v,m,b] = disc_dist(x,y)
    % Trim the spectrum for consistent x limits
    %idx1 = binary_search(x, -500);
    %idx2 = binary_search(x, 500);
    idx1 = 1;
    idx2 = size(x,1);
    
    x2 = x(idx1:idx2,1);
    ay = abs(y(idx1:idx2,1));
    ay = ay/sum(ay);
    
    v = sum(x2.*ay);
    b = sqrt(sum(x2.*x2.*ay) - v*v);
    m = max(abs(y));
end