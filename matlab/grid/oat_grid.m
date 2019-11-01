function [gvec, tc, gfuncs, sq2_base] = oat_grid(rootdir, n)

TOL = 1e-10;

% Assumes all files share the same time grid

% Get the base file
C = dlmread('/Users/tuckerkj/output/20190930/cumulant/vgs_o0/cumulant_N2000_D0p0_g2p_o0_f1p0_faa1.csv');
[tc, esc, essc] = unpack_symm(C);
sq2_base = get_squeezing(n, esc, essc);
%sq2_base_db = 10*log10(real(sq2_base));

if max(abs(imag(sq2_base))) > TOL
    disp('ERROR - Found nonzero imaginary part in sq2');
    return;
end

sq2_base = real(sq2_base);

filenames = dir([rootdir, '/*cumulant_N*.csv']);
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_gaa[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa[^_]+_gaa[^_]+.txt';
pat = '[^_]*cumulant_N[^_]+_D[^_]+_g[^_]+_o[^_]+_f[^_]+_faa([^_]+).csv';

A = zeros(size(filenames,1),2);
for i=1:size(filenames,1)
    t = regexp(filenames(i,1).name, pat, 'tokens');
    A(i,1) = str2double(strrep(t{1}(1), 'p', '.'));
    A(i,2) = i;
    filepaths{i} = [rootdir, '/', filenames(i,1).name];
end

A = sortrows(A,1);
gvec = unique(A(:,1));

gfuncs = zeros(size(tc,1),size(gvec,1));
for i=1:size(A,1)
    gidx = find(gvec == A(i,1));
    g = A(i,1);

    B = dlmread(filepaths{A(i,2)});
    %[es, ess] = unpack_symm_interp(B, tau);
    [tau, es, ess] = unpack_symm(B);
    sq2 = get_squeezing(n, es, ess);
    
    if max(abs(imag(sq2))) > TOL
        disp('ERROR - Found nonzero imaginary part in sq2');
        return;
    end
    
    sq2 = real(sq2);
    
    if norm(tc - tau) > TOL
        disp(['ERROR - Time grids do not align when gs = ', num2str(g)]);
        return;
    end
    
    gfuncs(:,gidx) = sq2 - sq2_base;
    %gfuncs(:,gidx) = 10*log10(sq2) - sq2_base_db;
    
    disp(['g = ', num2str(g)]);
end
