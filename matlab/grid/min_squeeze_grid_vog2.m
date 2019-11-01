function [gvec, ovec, minvals, finalvals, durs, sstimes] = min_squeeze_grid_vog2(rootdir, n, orvec, comp_ss)

filenames = dir([rootdir, '/*cumulant_N*.csv']);
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_gaa[^_]+.txt';
%pat = '[^_]*symm_N([^_]+)_D[^_]+_g[^_]+_o([^_]+)_f[^_]+_faa[^_]+_gaa[^_]+.txt';
pat = '[^_]*cumulant_N[^_]+_D[^_]+_g([^_]+)_o([^_]+)_f[^_]+_faa[^_]+.csv';

A = zeros(size(filenames,1),2);
for i=1:size(filenames,1)
    t = regexp(filenames(i,1).name, pat, 'tokens');
    A(i,1) = str2double(strrep(t{1}(1), 'p', '.'));
    A(i,2) = str2double(strrep(t{1}(2), 'p', '.'));
    A(i,3) = i;
    filepaths{i} = [rootdir, '/', filenames(i,1).name];
end

A = sortrows(A,1);
gvec = unique(A(:,1));
allovec = unique(A(:,2));

delo = 10;
mino = min(allovec);
maxo = max(allovec);
ovec = (mino:delo:maxo)';

% Hardcode f for now
f = 1;

% Initialize data by g value
data = cell(size(gvec,1),1);
for i=1:size(gvec,1)
    data{i,1}.ovals = zeros(size(orvec,1),1);
    data{i,1}.mvs = zeros(size(orvec,1),1);
end

% Populate data
for i=1:size(A,1)
    gidx = find(gvec == A(i,1));
    g = A(i,1);
    o = A(i,2);
    
    oc = (n/2)*sqrt(f^2 + g^2);
    or = round(o/oc,2);
    oidx = find(round(orvec,2) == or);
    
    if(size(oidx,1) == 0)
        if(size(orvec,1) == 1)
            oidx = 1;
        else
            disp('ERROR: Could not find Omega ratio in orvec')
            return
        end
    end
    
    data{gidx,1}.ovals(oidx,1) = o;

    B = dlmread(filepaths{A(i,3)});
    %[es, ess] = unpack_symm_interp(B, tau);
    [tau, es, ess] = unpack_symm(B);
    sq2 = get_squeezing(n, es, ess);
    if nargin > 3 && comp_ss ~= 0
        sq2ss = get_ss_sq2(n,o,f,g);
        [data{gidx,1}.mvs(oidx,1), finalvals(gidx,oidx,1), ~, durs(gidx,oidx,1), sstimes(gidx,oidx,1)] = get_min_sq2(tau, sq2, sq2ss);
    else
        %[data{gidx,1}.mvs(oidx,1), finalvals(gidx,oidx,1), ~, durs(gidx,oidx,1), sstimes(gidx,oidx,1)] = get_min_sq2(tau, sq2);
        data{gidx,1}.mvs(oidx,1) = get_min_sq2(tau, sq2);
    end
    
    disp(['g = ', num2str(g), ' o = ', num2str(orvec(oidx)), ' oidx = ', num2str(oidx)]);
end

% Interpolate data to fill in output structures one row at a time
minvals = zeros(size(gvec,1),size(ovec,1));
finalvals = zeros(size(gvec,1),size(ovec,1),1);
durs = zeros(size(gvec,1),size(ovec,1),1);
sstimes = zeros(size(gvec,1),size(ovec,1),1);
for gidx=1:size(gvec,1)
    for oidx=1:size(ovec,1)
        if size(data{gidx,1}.ovals,1) > 1
            o = ovec(oidx,1);
            if o >= min(data{gidx,1}.ovals) && o <= max(data{gidx,1}.ovals)
                minvals(gidx,oidx) = interp1(data{gidx,1}.ovals, data{gidx,1}.mvs, o);
            %else
            %    disp('ERROR: Found Omega value outside boundaries')
            end
        else
            minvals(gidx,oidx) = data{gidx,1}.mvs(1);
        end
    end
end

% minvals = zeros(size(gvec,1),size(orvec,1));
% finalvals = zeros(size(gvec,1),size(orvec,1),1);
% durs = zeros(size(gvec,1),size(orvec,1),1);
% sstimes = zeros(size(gvec,1),size(orvec,1),1);
% for gidx=1:size(gvec,1)
%     g = gvec(gidx);
%     oc = (n/2)*sqrt(f^2 + g^2);
%     for oidx=1:size(orvec,1)
%         or = orvec(oidx,1);
%         o = or*oc;
%         if o >= min(data{gidx,1}.ovals) && o <= max(data{gidx,1}.ovals)
%             minvals(gidx,oidx) = interp1(data{gidx,1}.ovals, data{gidx,1}.mvs, o);
%         end
%     end
% end