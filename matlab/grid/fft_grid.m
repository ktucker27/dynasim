function [V, M, B] = fft_grid(tt)

V = zeros(size(tt,1), size(tt,2));
M = zeros(size(tt,1), size(tt,2));
B = zeros(size(tt,1), size(tt,2));

for i=1:size(tt,1)
    for j=1:size(tt,2)
        [x,y] = plot_fft(tt{i,j}(:,1), tt{i,j}(:,2), 0);
        
        [V(i,j), M(i,j), B(i,j)] = disc_dist(x,y);
    end
end
end

function [v,m,b] = smooth_max(x,y)
    ay = smooth(abs(y), 11);
        
    [maxval, maxidx] = max(ay);
    
    v = x(maxidx);
    m = maxval;
    b = whm(x,ay);
end

function [v,m,b] = disc_dist(x,y)
    ay = abs(y);
    
    ay = ay/sum(ay);
    
    v = sum(x.*ay);
    b = sqrt(sum(x.*x.*ay) - v*v);
    m = max(smooth(abs(y), 11));
end