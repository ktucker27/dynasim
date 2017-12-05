function [V, M, B, T, G] = fft_grid(tt, tol)

if nargin == 1
    tol = 0.001;
end

V = zeros(size(tt,1), size(tt,2));
M = zeros(size(tt,1), size(tt,2));
B = zeros(size(tt,1), size(tt,2));
T = zeros(size(tt,1), size(tt,2));

for i=1:size(tt,1)
    for j=1:size(tt,2)
        [x,y] = plot_fft(tt{i,j}(:,1), tt{i,j}(:,2), 0);
        
        [V(i,j), M(i,j), B(i,j)] = disc_dist(x,y);
        T(i,j) = decay_time(tt{i,j}, tol);
    end
end

G = -log(tol)./T;

end

function [v,m,b] = smooth_max(x,y)
    ay = smooth(abs(y), 11);
        
    [maxval, maxidx] = max(ay);
    
    v = x(maxidx);
    m = maxval;
    b = whm(x,ay);
end

function [v,m,b] = disc_dist(x,y)
    % Trim the spectrum for consistent x limits
    idx1 = binary_search(x, -500);
    idx2 = binary_search(x, 500);
    
    x2 = x(idx1:idx2,1);
    ay = abs(y(idx1:idx2,1));
    ay = ay/sum(ay);
    
    v = sum(x2.*ay);
    b = sqrt(sum(x2.*x2.*ay) - v*v);
    m = max(abs(y));
end