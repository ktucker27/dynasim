function [mi, M, decay, freq, T] = exact_mi_grid(gvec, nvec, wvec, f, gamma, tau, savefile)

num_evals = 5;
mi = cell(size(nvec,1), size(gvec,1));
M = zeros(size(nvec,1), size(gvec,1));
T = zeros(size(nvec,1), size(gvec,1));
decay = zeros(size(nvec,1), size(gvec,1), num_evals);
freq = zeros(size(nvec,1), size(gvec,1), num_evals);

tol = 1.0e-12;

rho0 = get_rho0(nvec(1), pi/2, 0);

for i=1:size(gvec)
    disp(['g = ', num2str(gvec(i))]);
    for j=1:size(nvec)
        disp(['n = ', num2str(nvec(j))]);
        
        n = nvec(j);
        %rho0 = (1/(2^n))*ones(2^n,2^n);
        %rho0 = zeros(2^n,2^n);
        %rh0(1,1) = 1;
        %rho0 = get_rho0(n, pi/2);
        
        L = full(master_matrix(n, wvec(j), zeros(n,1), f, gvec(i), gamma));
        
        l = eig(L);
        lv = [real(l) imag(l)];
        lv = sortrows(lv, -1);
        num_found = 0;
        for lvidx = 1:size(lv,1)
            if abs(lv(lvidx, 2)) > tol
                decay(j,i,num_found+1) = abs(lv(lvidx,1));
                freq(j,i,num_found+1) = abs(lv(lvidx,2));
                if num_found == num_evals
                    break;
                end
                num_found = num_found + 1;
            end
        end
        
        if lv(1,1) > tol || lv(1,2) > tol
            disp('WARNING: Found nonzero first eigenvalue of L');
        end
        
        % Since g is to be plotted on the x-axis of a color plot, it varies
        % with the columns of the return matrices
        mi{j,i} = [tau, mutual_info_vt(rho0,L,tau)];
        
        [M(j,i), maxidx] = max(abs(mi{j,i}(:,2)));
        T(j,i) = mi{j,i}(maxidx,1);
        
        if nargin > 7
            save(savefile, 'gvec', 'nvec', 'wvec', 'f', 'gamma', 'rho0', 'mi', 'M', 'T');
        end
    end
end