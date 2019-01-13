function [ovec, minvals, finalvals, durs] = min_squeeze_grid_mcwf(rootdir, n)

filenames = dir([rootdir, '/o*']);
pat = 'o([^_]+)';

A = zeros(size(filenames,1),2);
for i=1:size(filenames,1)
    t = regexp(filenames(i,1).name, pat, 'tokens');
    A(i,1) = str2double(strrep(t{1}(1), 'p', '.'));
    A(i,2) = i;
    filepaths{i} = [rootdir, '/', filenames(i,1).name];
end

A = sortrows(A,1);
ovec = unique(A(:,1));

minvals = zeros(size(ovec,1),1);
finalvals = zeros(size(ovec,1),1);
durs = zeros(size(ovec,1),1);
for i=1:size(A,1)
    oidx = find(ovec == A(i,1));
    
    [tau, es, ess] = proc_mcwf_output(filepaths{A(i,2)}, 0, 1);
    sq2 = get_squeezing(n, es, ess);
    [minvals(oidx,1), finalvals(oidx,1), ~, durs(oidx,1)] = get_min_sq2(tau, sq2);
    
    disp([' o = ', num2str(ovec(oidx))]);
end