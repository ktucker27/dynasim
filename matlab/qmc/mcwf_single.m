function [tau, rho, cdfs, outcomes, szs, es, ess, js, ms] = mcwf_single(n, w, o, faa, fab, chi, gamma, gel, c, t0, dt, t, num_trajectories)

dim = (n/2 + 1)^2;

gc = gamma*fab;
gl = gamma*(faa - fab);
kl = w;         % ? Double check in master equation
dl = gamma*gel; % ?
oms = chi;

dur = t - t0;
num_times = ceil(dur/dt);
tau = zeros(num_times,1);
rho = zeros(dim, dim, num_times);
cdfs = zeros(4,num_times-1);
outcomes = zeros(num_times-1,1);
es = zeros(3, num_times);
ess = zeros(3, 3, num_times);

tau(1,1) = t0;

szs = zeros(num_times,1);

% TODO - Assumes spin up initial conditions
szs(1,1) = n/2;
js = zeros(num_trajectories,num_times);
ms = zeros(num_trajectories,num_times);
c = zeros(num_trajectories, 1);
for i = 1:num_trajectories
    c(i,1) = 1;
    js(i,1) = n/2;
    ms(i,1) = n/2;
end

% Iterate over time points
for ti = 2:num_times
    tau(ti,1) = t0 + (ti-1)*dt;
    
    if mod(ti-1,100) == 0
        disp(['t = ', num2str(tau(ti,1))]);
    end
    
    % Iterate over the trajectories updating each one and computing the
    % expected values for this time point
    for i = 1:num_trajectories
        % Update the trajectory
        Jt = js(i,ti - 1);
        Mt = ms(i,ti - 1);
        
        % Determine the jump
        cdf = get_cdf(c(i), n, Jt, Mt, gc, gl, dt);
        outcome = roll_dice(cdf);
                
        cdfs(:,ti-1,i) = cdf;
        outcomes(ti-1,i) = outcome;
        
        if outcome == 1
            % Collective decay
            [~, ajmm] = get_ajm(Jt,Mt);
            c(i,1) = sqrt(dt*gc)*ajmm*c(i,1);
            Mt = Mt - 1;
        elseif outcome == 2
            % Individual decay, s = 0
            [~, ajmm] = get_ajm(Jt, Mt);
            [p0, ~, ~] = get_pjms(n, Jt, Mt, ajmm);
            c(i,1) = sqrt(dt*gl)*p0*c(i,1);
            Mt = Mt - 1;
        elseif outcome == 3
            % Individual decay, s = -1
            s = -1;
            
            [~, ajmm] = get_ajm(Jt, Mt);
            [~, pm, ~] = get_pjms(n, Jt, Mt, ajmm);
            c(i,1) = sqrt(dt*gl)*pm*c(i,1);
            
            Jt = Jt + s;
            Mt = Mt - 1;
        elseif outcome == 4
            % Individual decay, s = 1
            s = 1;
            
            [~, ajmm] = get_ajm(Jt, Mt);
            [~, ~, pp] = get_pjms(n, Jt, Mt, ajmm);
            c(i,1) = sqrt(dt*gl)*pp*c(i,1);
            
            Jt = Jt + s;
            Mt = Mt - 1;
        else
            % No jump
            [~, ajmm] = get_ajm(Jt, Mt);
            omjm = oms*ajmm*ajmm ...
                - 1i*0.5*(gc*ajmm*ajmm + gl*(n/2 + Mt) + kl*(n/2 - Mt) + dl*n);
            c(i,1) = (1 - 1i*omjm*dt)*c(i,1);
        end
        c(i,1) = c(i,1)/abs(c(i,1));
        
        js(i,ti) = Jt;
        ms(i,ti) = Mt;
        
        szs(ti,1) = szs(ti,1) + Mt;
    end
    
    % Normalize the density
    szs(ti,1) = real(szs(ti,1)/num_trajectories);
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
if J == 0 && M == 0
    p0 = 0;
    pm = 0;
else
    p0 = sqrt((2+n)/(4*J*(J+1)))*ajmm;
    pm = -sqrt((n + 2*J + 2)*(J + M)*(J + M - 1)/(4*J*(2*J + 1)));
end
pp = sqrt((n - 2*J)*(J - M + 1)*(J - M + 2)/(4*(J + 1)*(2*J + 1)));
end

% Returns the jump probability CDF for use in roll_dice
function cdf = get_cdf(cjm, n, J, M, gc, gl, dt)
cdf = zeros(4,1);

dnj = get_dnj(n,J);
[~, ajmm] = get_ajm(J,M);

% Collective decay
cdf(1,1) = cdf(1,1) + abs(ajmm*cjm)^2;

% Individual decay
[p0, pm, pp] = get_pjms(n, J, M, ajmm);
cdf(2,1) = cdf(2,1) + abs(p0*cjm)^2;
cdf(3,1) = cdf(3,1) + abs(pm*cjm)^2;
cdf(4,1) = cdf(4,1) + abs(pp*cjm)^2;

cdf(1,1) = dt*gc*dnj*cdf(1,1);
cdf(2,1) = dt*gl*dnj*cdf(2,1);
cdf(3,1) = dt*gl*dnj*cdf(3,1);
cdf(4,1) = dt*gl*dnj*cdf(4,1);

cdf = cdf/dnj;

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