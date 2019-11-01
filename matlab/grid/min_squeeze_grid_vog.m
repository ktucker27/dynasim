function [gvec, ovec, minvals, finalvals, durs, sstimes, mintimes] = min_squeeze_grid_vog(rootdir, n, comp_ss)

filenames = dir([rootdir, '/*cumulant_N*.csv']);
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_gaa[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa[^_]+_gaa[^_]+.txt';
pat = '[^_]*cumulant_N[^_]+_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa([^_]+).csv';

A = zeros(size(filenames,1),2);
for i=1:size(filenames,1)
    t = regexp(filenames(i,1).name, pat, 'tokens');
    A(i,2) = str2double(strrep(t{1}(1), 'p', '.'));
    A(i,1) = str2double(strrep(t{1}(2), 'p', '.'));
    A(i,3) = i;
    filepaths{i} = [rootdir, '/', filenames(i,1).name];
end

A = sortrows(A,1);
gvec = unique(A(:,1));
ovec = unique(A(:,2));

% Hardcode f for now
f = 1;

minvals = zeros(size(gvec,1),size(ovec,1),1);
finalvals = zeros(size(gvec,1),size(ovec,1),1);
durs = zeros(size(gvec,1),size(ovec,1),1);
sstimes = zeros(size(gvec,1),size(ovec,1),1);
mintimes = zeros(size(gvec,1),size(ovec,1),1);
for i=1:size(A,1)
    gidx = find(gvec == A(i,1));
    oidx = find(ovec == A(i,2));
    g = A(i,1);
    o = A(i,2);

    B = dlmread(filepaths{A(i,3)});
    %[es, ess] = unpack_symm_interp(B, tau);
    [tau, es, ess] = unpack_symm(B);
    sq2 = get_squeezing(n, es, ess);
    if nargin > 2 && comp_ss ~= 0
        sq2ss = get_ss_sq2(n,o,f,g);
        [minvals(gidx,oidx,1), finalvals(gidx,oidx,1), mintimes(gidx, oidx, 1), durs(gidx,oidx,1), sstimes(gidx,oidx,1)] = get_min_sq2(tau, sq2, sq2ss);
    else
        [minvals(gidx,oidx,1), finalvals(gidx,oidx,1), mintimes(gidx, oidx, 1), durs(gidx,oidx,1), sstimes(gidx,oidx,1)] = get_min_sq2(tau, sq2, 0, 0);
    end
    
    disp(['g = ', num2str(g), ' o = ', num2str(ovec(oidx))]);
end