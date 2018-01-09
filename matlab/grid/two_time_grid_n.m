function [gvec,nvec,tt] = two_time_grid_n(rootdir)

filenames = dir([rootdir, '/*time_corr*.txt']);
pat = '[^_]*time_corr_N([^_]+)_D[^_]+_g([^_]+)_w([^_]+)_f[^_]+.txt';

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

tt = cell(size(nvec,1), size(gvec,1));
for i=1:size(A,1)
    nidx = find(nvec == A(i,1));
    gidx = find(gvec == A(i,2));

    f = dlmread(filepaths{A(i,3)});
    tt{nidx, gidx} = [f(:,1),f(:,2) + 1i*f(:,3)];
    
    disp(['g = ', num2str(gvec(gidx)), ', n = ', num2str(nvec(nidx))]);
end