function A = coll_two_time_matrix_to(n, w, f, g, gamma)

cz = coll_steady_vals(n, w, f, gamma);

A(1,1) = -0.5*(gamma + w) + 0.5*(f - 1i*g)*(n-2)*cz;
A(1,2) = 0.5*(f - 1i*g)*cz;

A(2,1) = 0.5*(f - 1i*g)*(n - 1)*cz;
A(2,2) = -0.5*(gamma + w);