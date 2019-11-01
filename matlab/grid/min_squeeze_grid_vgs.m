function [gsvec, ovec, minvals, finalvals, durs, sstimes, ssvals, tfs] = min_squeeze_grid_vgs(rootdir, n, f, g)

filenames = dir([rootdir, '/*cumulant_N*.csv']);
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_gaa[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa[^_]+_gaa[^_]+.txt';
pat = '[^_]*cumulant_N[^_]+_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa([^_]+).csv';

A = zeros(size(filenames,1),2);
for i=1:size(filenames,1)
    t = regexp(filenames(i,1).name, pat, 'tokens');
    A(i,1) = str2double(strrep(t{1}(1), 'p', '.'));
    A(i,2) = str2double(strrep(t{1}(2), 'p', '.'));
    A(i,3) = i;
    filepaths{i} = [rootdir, '/', filenames(i,1).name];
end

A = sortrows(A,2);
ovec = unique(A(:,1));
gsvec = unique(A(:,2)) - 1;

minvals = zeros(size(ovec,1),1);
finalvals = zeros(size(ovec,1),1);
durs = zeros(size(ovec,1),1);
sstimes = zeros(size(ovec,1),1);
ssvals = zeros(size(ovec,1),3);
tfs = zeros(size(ovec,1),1);
for i=1:size(A,1)
    o = A(i,1);
    gs = A(i,2) - 1;
    gsidx = find(gsvec == gs);
    gsvec(gsidx,1) = gs;

    tmin = min([0.05*10^(-floor(log10(gs))), .01]);

    B = dlmread(filepaths{A(i,3)});
    %[es, ess] = unpack_symm_interp(B, tau);
    [tau, es, ess] = unpack_symm(B);
    sq2 = get_squeezing(n, es, ess);
    ssvals(gsidx,:) = es(:,end)';
    tfs(gsidx,1) = tau(end);
    if nargin > 2
        sq2ss = get_ss_sq2(n,o,f,g);
        [minvals(gsidx,1), finalvals(gsidx,1), ~, durs(gsidx,1), sstimes(gsidx,1)] = get_min_sq2(tau, sq2, sq2ss, tmin);
    else
        [minvals(gsidx,1), finalvals(gsidx,1), ~, durs(gsidx,1), sstimes(gsidx,1)] = get_min_sq2(tau, sq2, 0, tmin);
    end
    
    disp(['o = ', num2str(o), ' gs = ', num2str(gs)]);
end