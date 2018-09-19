function [tau, rho, cdfs, outcomes, szs, es, ess] = mcwf(n, w, o, faa, fab, chi, gamma, gel, c, t0, dt, t, num_trajectories)

dim = size(c,1);

[sp, sm, sz] = get_ops(n);

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

% Update the density operator
tau(1,1) = t0;

idx1 = 1;
for ji = n/2:-1:-n/2
    jstart = idx1;
    for M1 = ji:-1:-ji
        idx2 = jstart;
        for M2 = ji:-1:-ji
            rho(idx1, idx2, 1) = rho(idx1, idx2, 1) + c(idx1,1)*c(idx2,1)';
            idx2 = idx2 + 1;
        end
        idx1 = idx1 + 1;
    end
end

szs = zeros(num_times,1);
szs(1,1) = trace(rho(:,:,1)*sz);

[es, ess] = update_evs(c, sp, sm, sz, es, ess, 1);

J = zeros(num_trajectories,1);
cs = zeros(dim, num_trajectories);
for i = 1:num_trajectories
    cs(:,i) = c;
    J(i,1) = n/2;
end

% Iterate over time points
for ti = 2:num_times
    tau(ti,1) = t0 + (ti-1)*dt;
    
    if mod(ti-1,100) == 0
        disp(['t = ', num2str(tau(ti,1))]);
    end
    
    % Iterate over the trajectories updating each one and computing the
    % density matrix for this time point
    for i = 1:num_trajectories
        % Update the trajectory
        Jt = J(i,1);
        
        % Determine the jump
        cdf = get_cdf(cs(:,i), n, Jt, gc, gl, dt);
        outcome = roll_dice(cdf);
                
        if i == 1
            cdfs(:,ti-1) = cdf;
            outcomes(ti-1,1) = outcome;
        end
        
        % Update the state coefficients
        newc = zeros(dim,1);
        
        idx = 1;
        for ji = n/2:-1:Jt+1
            idx = idx + 2*ji + 1;
        end
        
        if outcome == 1
            % Collective decay
            for M = Jt:-1:-Jt + 1
                [~, ajmm] = get_ajm(Jt,M);
                newc(idx + 1,1) = sqrt(dt*gc)*ajmm*cs(idx,i);
                idx = idx + 1;
            end
        elseif outcome == 2
            % Individual decay, s = 0
            for M = Jt:-1:-Jt + 1
                [~, ajmm] = get_ajm(Jt, M);
                [p0, ~, ~] = get_pjms(n, Jt, M, ajmm);
                newc(idx + 1,1) = sqrt(dt*gl)*p0*cs(idx,i);
                idx = idx + 1;
            end
        elseif outcome == 3
            % Individual decay, s = -1
            s = -1;
            newidx = 1;
            for ji = n/2:-1:(Jt + s) + 1
                newidx = newidx + 2*ji + 1;
            end
            
            for M = Jt:-1:-Jt+2
                [~, ajmm] = get_ajm(Jt, M);
                [~, pm, ~] = get_pjms(n, Jt, M, ajmm);
                newc(newidx,1) = sqrt(dt*gl)*pm*cs(idx,i);
                idx = idx + 1;
                newidx = newidx + 1;
            end
            
            Jt = Jt + s;
        elseif outcome == 4
            % Individual decay, s = 1
            s = 1;
            newidx = 1;
            for ji = n/2:-1:(Jt + s) + 1
                newidx = newidx + 2*ji + 1;
            end
            
            % Advance the new index to start one less than the original M
            newidx = newidx + 2;
            
            for M = Jt:-1:-Jt
                [~, ajmm] = get_ajm(Jt, M);
                [~, ~, pp] = get_pjms(n, Jt, M, ajmm);
                newc(newidx,1) = sqrt(dt*gl)*pp*cs(idx,i);
                idx = idx + 1;
                newidx = newidx + 1;
            end
            
            Jt = Jt + s;
        else
            % No jump
            for M = Jt:-1:-Jt
                [~, ajmm] = get_ajm(Jt, M);
                omjm = oms*ajmm*ajmm ...
                     - 1i*0.5*(gc*ajmm*ajmm + gl*(n/2 + M) + kl*(n/2 - M) + dl*n);
                newc(idx,1) = (1 - 1i*omjm*dt)*cs(idx,i);
                idx = idx + 1;
            end
        end
        cs(:,i) = newc/norm(newc);
        
        % Update the density operator
%         rho_traj = zeros(dim,dim);
%         idx1 = 1;
%         for ji = n/2:-1:-n/2
%             jstart = idx1;
%             for M1 = ji:-1:-ji
%                 idx2 = jstart;
%                 for M2 = ji:-1:-ji
%                     rho_traj(idx1, idx2) = cs(idx1,i)*cs(idx2,i)';
%                     rho(idx1, idx2, ti) = rho(idx1, idx2, ti) + cs(idx1,i)*cs(idx2,i)';
%                     idx2 = idx2 + 1;
%                 end
%                 idx1 = idx1 + 1;
%             end
%         end
        
        %szs(ti,1) = szs(ti,1) + trace(rho_traj*sz);
        szs(ti,1) = szs(ti,1) + cs(:,i)'*sz*cs(:,i);
        [es, ess] = update_evs(cs(:,i), sp, sm, sz, es, ess, ti);
        J(i,1) = Jt;
    end
    
    % Normalize the density
    %rho(:,:,ti) = rho(:,:,ti)/num_trajectories;
    szs(ti,1) = real(szs(ti,1)/num_trajectories);
    es(:,ti) = es(:,ti)/num_trajectories;
    ess(:,:,ti) = ess(:,:,ti)/num_trajectories;
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
    cdf(3,1) = cdf(3,1) + abs(pm*cjm)^2;
    cdf(4,1) = cdf(4,1) + abs(pp*cjm)^2;

    idx = idx + 1;
end

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

function [sp, sm, sz] = get_ops(n)
dim = (n/2 + 1)^2;
sp = zeros(dim, dim);
sm = zeros(dim, dim);
sz = zeros(dim, dim);
idx = 1;
for ji = n/2:-1:0
    mvec = ji:-1:-ji;
    ms = size(mvec,2);
    sz(idx:idx + ms - 1, idx:idx + ms - 1) = diag(mvec);
    sp(idx:idx + ms - 1, idx:idx + ms - 1) = circshift(sqrt(diag((n/2 - mvec).*(n/2 + mvec + 1))),-1);
    sm(idx:idx + ms - 1, idx:idx + ms - 1) = circshift(sqrt(diag((n/2 + mvec).*(n/2 - mvec + 1))),1);
    idx = idx + ms;
end
end

function [es, ess] = update_evs(c, sp, sm, sz, es, ess, i)
sx = 0.5*(sp + sm);
sy = 1i*0.5*(sm - sp);

es(1,i) = es(1,i) + c'*sx*c;
es(2,i) = es(2,i) + c'*sy*c;
es(3,i) = es(3,i) + c'*sz*c;
end