function [tau, rho] = mcwf(n, w, o, faa, fab, chi, gamma, gel, c, t0, dt, t, num_trajectories)

J = 1;

dim = size(c,1);

gc = gamma*fab;
gl = gamma*(faa - fab);
kl = w;         % ? Double check in master equation
dl = gamma*gel; % ?
oms = chi;

dur = t - t0;
num_times = ceil(dur/dt);
tau = zeros(num_times,1);
rho = zeros(dim, dim, num_times);

% Iterate over time points
for ti = 1:num_times
    tau(ti,1) = t0 + ti*dt;
    
    if mod(ti,100) == 0
        disp(['t = ', num2str(tau(ti,1))]);
    end
    
    % Iterate over the trajectories updating each one and computing the
    % density matrix for this time point
    for i = 1:num_trajectories
        % Update the trajectory
        
        % Determine the jump
        cdf = get_cdf(c, J, gc, gl, dt);
        outcome = roll_dice(cdf);
        
        % Update the state coefficients
        newc = zeros(dim,1);
        
        idx = 1;
        for ji = n/2:-1:J+1
            idx = idx + 2*ji + 1;
        end
        
        if outcome == 1
            % Collective decay
            for M = J:-1:-J + 1
                [~, ajmm] = get_ajm(J,M);
                newc(idx + 1,1) = sqrt(dt*gc)*ajmm*c(idx,1);
                idx = idx + 1;
            end
        elseif outcome == 2
            % Individual decay, s = 0
            for M = J:-1:-J
                [~, ajmm] = get_ajm(J, M);
                [p0, ~, ~] = get_pjms(n, J, M, ajmm);
                newc(idx + 1,1) = sqrt(dt*gl)*p0*c(idx,1);
                idx = idx + 1;
            end
        elseif outcome == 3
            % Individual decay, s = -1
            s = -1;
            newidx = 1;
            for ji = n/2:-1:(J + s) + 1
                newidx = newidx + 2*ji + 1;
            end
            
            for M = J:-1:-J+1
                [~, ajmm] = get_ajm(J, M);
                [~, pm, ~] = get_pjms(n, J, M, ajmm);
                newc(newidx,1) = sqrt(dt*gl)*pm*c(idx,1);
                idx = idx + 1;
                newidx = newidx + 1;
            end
            
            J = J + s;
        elseif outcome == 4
            % Individual decay, s = 1
            s = 1;
            newidx = 1;
            for ji = n/2:-1:(J + s) + 1
                newidx = newidx + 2*ji + 1;
            end
            
            % Advance the new index to start one less than the original M
            newidx = newidx + 2;
            
            for M = J:-1:-J
                [~, ajmm] = get_ajm(J, M);
                [~, ~, pp] = get_pjms(n, J, M, ajmm);
                newc(newidx,1) = sqrt(dt*gl)*pp*c(idx,1);
                idx = idx + 1;
                newidx = newidx + 1;
            end
            
            J = J + s;
        else
            % No jump
            for M = J:-1:-J
                [~, ajmm] = get_ajm(J, M);
                omjm = oms*ajmm*ajmm ...
                     - 1i*0.5*(gc*ajmm*ajmm + gl*(n/2 + M) + kl*(n/2 - M) + dl*n);
                newc(idx,1) = (1 - 1i*omjm*dt)*c(idx,1);
                idx = idx + 1;
            end
        end
        c = newc;
        
        % Update the density operator
        idx1 = 1;
        for J = n/2:-1:-n/2
            jstart = idx1;
            for M1 = J:-1:-J
                idx2 = jstart;
                for M2 = J:-1:-J
                    rho(idx1, idx2, ti) = rho(idx1, idx2, ti) + c(idx1,1)*c(idx2,1)';
                    idx2 = idx2 + 1;
                end
                idx1 = idx1 + 1;
            end
        end
    end
    
    % Normalize the density
    rho(:,:,ti) = rho(:,:,ti)/num_trajectories;
end
end

function dnj = get_dnj(n,J)
dnj = factorial(n)*(2*J+1)/(factorial(n/2 - J)*factorial(n/2 + J + 1));
end

function [ajmp, ajmm] = get_ajm(J,M)
ajmp = sqrt((J-M)*(J+M+1));
ajmm = sqrt((J+M)*(J-M+1));
end

function [p0, pm, pp] = get_pjms(n, J, M, ajmm)
p0 = sqrt((2+n)/(4*J*(J+1)))*ajmm;
pm = -sqrt((n + 2*J + 2)*(J + M)*(J + M - 1)/(4*J*(2*J + 1)));
pp = sqrt((n - 2*J)*(J - M + 1)*(J - M + 2)/(4*(J + 1)*(2*J + 1)));
end

% Returns the jump probability CDF for use in roll_dice
function cdf = get_cdf(c, n, J, gc, gl, dt)
cdf = zeros(4,1);

dnj = get_dnj(n,J);

% Advance the index to the start of the Dicke manifold corresponding to J
idx = 1;
for ji = n/2:-1:J+1
    idx = idx + 2*ji + 1;
end

for M = J:-1:-J
    cjm = c(idx,1);
    [~, ajmm] = get_ajm(J,M);
    
    % Collective decay
    cdf(1,1) = cdf(1,1) + abs(ajmm*cjm)^2;
    
    % Individual decay
    [p0, pm, pp] = get_pjms(n, J, M, ajmm);
    cdf(2,1) = cdf(2,1) + abs(p0*cjm)^2;
    cdf(3,1) = cdf(2,1) + abs(pm*cjm)^2;
    cdf(4,1) = cdf(2,1) + abs(pp*cjm)^2;

    idx = idx + 1;
end

cdf(1,1) = dt*gc*dnj*cdf(1,1);
cdf(2,1) = dt*gl*dnj*cdf(1,1);
cdf(3,1) = dt*gl*dnj*cdf(1,1);
cdf(4,1) = dt*gl*dnj*cdf(1,1);

% Go from probabilities to CDF
for i=2:size(cdf,1)
    cdf(i,1) = cdf(i-1,1) + cdf(i,1);
end

end

function idx = roll_dice(cdf)
eps = rand();
idx = 1;
for i=1:size(cdf,1)
    if eps <= cdf(i,1)
        break;
    end
    
    idx = idx + 1;
end
end