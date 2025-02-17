function A = coll_two_time_matrix(n, w, f, g, gamma)
% -------------------------------------------------------------------------
% Returns the matrix governing the collective two-time correlation
% equations using a fourth order cumulant expansion.  The ordering of the
% variables is as follows:
%
% <\sigma_a^+(t + tau) sigma_b^-(t)>
% <\sigma_a^+(t + tau) sigma_a^-(t)>
% <\sigma_a^z(t + tau) sigma_b^+(t + tau) sigma_c^-(t)>
% <\sigma_a^z(t + tau) sigma_b^+(t + tau) sigma_b^-(t)>
% <\sigma_a^z(t + tau) sigma_b^+(t + tau) sigma_a^-(t)>
% -------------------------------------------------------------------------

[cz, czz, cpm] = coll_steady_vals(n, w, f, gamma);

A = zeros(5,5);

A(1,1) = -0.5*(gamma + w);
A(1,3) = 0.5*gamma*(f - 1i*g)*(n - 2);
A(1,4) = 0.5*gamma*(f - 1i*g);

A(2,2) = -0.5*(gamma + w);
A(2,5) = 0.5*gamma*(f - 1i*g)*(n - 1);

A(3,1) = -(gamma - w) - 0.5*gamma*(f + 1i*g) ...
    + 0.5*gamma*(f - 1i*g)*(n - 3)*(czz - 2*cz*cz) ...
    -2*gamma*(f + 1i*g)*(n - 2)*cpm ...
    -2*gamma*(f - 1i*g)*(n - 5/2)*cpm;
A(3,2) = 0.5*gamma*(f - 1i*g)*(czz - 2*cz*cz) ...
    - gamma*(f - 1i*g)*cpm;
A(3,3) = -1.5*(gamma + w) - gamma*f ...
    + gamma*(f - 1i*g)*(n - 3)*cz;
A(3,4) = gamma*(f - 1i*g)*cz;

A(4,1) = -0.5*gamma*(f + 1i*g) ...
    + 0.5*gamma*(f - 1i*g)*(n - 2)*(czz - 2*cz*cz) ...
    - 2*gamma*f*(n - 2)*cpm;
A(4,2) = -(gamma - w) - 2*gamma*f*(n - 2)*cpm;
A(4,3) = 0.5*gamma*(f - 1i*g)*(n - 2)*cz;
A(4,4) = -1.5*(gamma + w);
A(4,5) = -gamma*f ...
    + 0.5*gamma*(f - 1i*g)*(n - 2)*cz;

A(5,1) = -(gamma - w) ...
    + 0.5*gamma*(f - 1i*g)*(n - 2)*(czz - 2*cz*cz) ...
    - gamma*(3*f - 1i*g)*(n - 2)*cpm;
A(5,2) = -0.5*gamma*(f + 1i*g) ...
    - gamma*(f + 1i*g)*(n - 2)*cpm;
A(5,3) = 0.5*gamma*(f - 1i*g)*(n - 2)*cz;
A(5,4) = -gamma*f;
A(5,5) = -1.5*(gamma + w) ...
    + 0.5*gamma*(f - 1i*g)*(n - 2)*cz;
