function rho0 = get_rho0(n, theta, phi)

rand_phase = (nargin < 3);

for i=1:n
    if rand_phase
        phi = (pi/6)*rand();
    end
    
    v = [cos(theta/2);exp(1i*phi)*sin(theta/2)];
    rhon = v*v';
    if i == 1
        rho0 = rhon;
    else
        rho0 = kron(rho0, rhon);
    end
end