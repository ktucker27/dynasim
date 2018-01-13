function [V, G] = coll_two_time_rates_grid(n, w, f, g, gamma)

V = zeros(size(n,1), size(g,1));
G = zeros(size(n,1), size(g,1));
for i=1:size(g,1)
    [G(:,i), V(:,i)] = coll_two_time_rates(n, w, f, g(i), gamma);
end