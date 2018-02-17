function [gvec,nvec,mi,M,T] = mi_grid(rootdir)

filenames = dir([rootdir, '/*exact_N*.txt']);
pat = '[^_]*exact_N([^_]+)_D[^_]+_g([^_]+)_w([^_]+)_f[^_]+.txt';

A = zeros(size(filenames,1),3);
for i=1:size(filenames,1)
    t = regexp(filenames(i,1).name, pat, 'tokens');
    A(i,1) = str2double(strrep(t{1}(1), 'p', '.'));
    A(i,2) = str2double(strrep(t{1}(2), 'p', '.'));
    A(i,3) = i;
    filepaths{i} = [rootdir, '/', filenames(i,1).name];
end

A = sortrows(A,[1,2]);
nvec = unique(A(:,1));
gvec = unique(A(:,2));

mi = cell(size(nvec,1), size(gvec,1));
M = zeros(size(nvec,1), size(gvec,1));
T = zeros(size(nvec,1), size(gvec,1));
for i=1:size(A,1)
    nidx = find(nvec == A(i,1));
    gidx = find(gvec == A(i,2));

    f = dlmread(filepaths{A(i,3)});
    pms = f(:,2);
    zs = f(:,5);
    zzs = f(:,8);
    
    ps = f(:,3) + 1i*f(:,4);
    zps = f(:,6) + 1i*f(:,7);
    pps = f(:,9) + 1i*f(:,10);
    
    mi{nidx, gidx} = mutual_info2(pms, zs, zzs, ps, zps, pps);
    [M(nidx, gidx), maxidx] = max(abs(mi{nidx, gidx}));
    T(nidx, gidx) = f(maxidx,1);
    
    disp(['g = ', num2str(gvec(gidx)), ', n = ', num2str(nvec(nidx))]);
end