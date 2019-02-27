function [nvec, ovec, minvals, finalvals, durs, sstimes] = min_squeeze_grid(rootdir, tau, f, g)

filenames = dir([rootdir, '/*cumulant_N*.csv']);
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_gaa[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa[^_]+_gaa[^_]+.txt';
pat = '[^_]*cumulant_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa[^_]+.csv';

A = zeros(size(filenames,1),2);
for i=1:size(filenames,1)
    t = regexp(filenames(i,1).name, pat, 'tokens');
    A(i,1) = str2double(strrep(t{1}(1), 'p', '.'));
    A(i,2) = str2double(strrep(t{1}(2), 'p', '.'));
    A(i,3) = i;
    filepaths{i} = [rootdir, '/', filenames(i,1).name];
end

A = sortrows(A,2);
ovec = unique(A(:,2));

nvec = zeros(size(ovec,1),1);
minvals = zeros(size(ovec,1),1);
finalvals = zeros(size(ovec,1),1);
durs = zeros(size(ovec,1),1);
sstimes = zeros(size(ovec,1),1);
for i=1:size(A,1)
    oidx = find(ovec == A(i,2));
    n = A(i,1);
    o = A(i,2);
    nvec(oidx,1) = n;
    
%     if n == 100
%         tau = (0:.001:5)';
%     end

    B = dlmread(filepaths{A(i,3)});
    %[es, ess] = unpack_symm_interp(B, tau);
    [tau, es, ess] = unpack_symm(B);
    sq2 = get_squeezing(n, es, ess);
    if nargin > 2
        sq2ss = get_ss_sq2(n,o,f,g);
        [minvals(oidx,1), finalvals(oidx,1), ~, durs(oidx,1), sstimes(oidx,1)] = get_min_sq2(tau, sq2, sq2ss);
    else
        [minvals(oidx,1), finalvals(oidx,1), ~, durs(oidx,1), sstimes(oidx,1)] = get_min_sq2(tau, sq2);
    end
    
    disp(['n = ', num2str(n), ' o = ', num2str(ovec(oidx))]);
end