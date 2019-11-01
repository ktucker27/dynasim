function [t, es, ess, sq2, j2, num_traj] = proc_mcwf_traj(rootdir, quiet)

% proc_mcwf_traj: Organizes results coming out of single trajectories from
%                 the MCWF Java executable. Assumes one trajectory per file
%
% INPUT: rootdir  - Can be a directory of files to aggregate, or a
%                   single file to process
%
% OUTPUT: t - Vector of times
%         es - Expected values of collective spins. es(:,i,j) is the
%              collective Bloch vector at time index i for trajectory j
%         ess - Second order expected values. ess(i,j,k,l) is the expected
%               value <J^i J^j> at time index k of trajectory l
%         sq2 - sq2(:,i) is the squeezing function for trajectory i
%         num_traj - Total number of trajectories found in the rootdir
%

if nargin < 2
    quiet = 0;
end

if ~isfolder(rootdir)
    b = split(rootdir, '/');
    files.name = b{size(b,1),1};
    c = join(b(1:end-1,1),'/');
    files.folder = c{1,1};
else
    files = dir([rootdir, '/*mcwf*.txt']);
end

num_traj = 0;
for i=1:size(files)
    filepath = [files(i,1).folder, '/', files(i,1).name];
    
    if ~quiet
        disp(['Processing: ', files(i,1).name]);
    end
    
    % Get header information
    header = dlmread(filepath, ',', [0 0 0 0]);
    if header ~= 1
        disp('ERROR - Found file with more than one trajectory');
        return;
    end
    num_traj = num_traj + header(1,1);
    
    % Process the data
    % TODO - Is there a way to do this with only one file read?
    data = dlmread(filepath, ',', 1, 0);
    
    [~, es_traj, ess_traj, ~, ~] = unpack_symm(data);
    
    % Grab or validate the time vector
    if i == 1
        t = data(:,1);
        n = round(2*sqrt(sum(es_traj(:,1).^2)));
        
        es = zeros(3, size(t,1), size(files,1));
        ess = zeros(3, 3, size(t,1), size(files,1));
        sq2 = zeros(size(t,1), size(files,1));
        j2 = zeros(size(t,1), size(files,1));
    else
        if norm(t - data(:,1)) ~= 0
            disp(['ERROR - Times in file ', files(i,1).name, ' do not align']);
            return;
        end
        
        if abs(n - round(2*sqrt(sum(es_traj(:,1).^2)))) ~= 0
            disp(['ERROR - Found value of n in file ', files(i,1).name, ' that does not match']);
            return;
        end
    end
    
    es(:,:,i) = es_traj;
    ess(:,:,:,i) = ess_traj;
    [sq2(:,i), ~, ~, ~] = get_squeezing(n, es_traj, ess_traj);
    
    j2(:,i) = get_j2(ess_traj);
end