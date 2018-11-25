function [t, es, ess, num_traj, esSe, essSe] = proc_mcwf_output(rootdir)

% proc_mcwf_output: Aggregates and organizes results coming out of
%                   the MCWF Java executable
%
% INPUT: rootdir - Can be a directory of files to aggregate, or a
%                  single file to process
%
% OUTPUT: t - Vector of times
%         es - Expected values of collective spins. es(:,i) is the
%              collective Bloch vector at time index i
%         ess - Second order expected values. ess(i,j,k) is the expected
%               value <J^i J^j> at time index k
%         num_traj - Total number of trajectories found in the aggregated
%                    files
%         esSe - Statistical error of the first order expectations
%         essSe - Statistical error of the second order expectations

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
    
    disp(['Processing: ', files(i,1).name]);
    
    % Get header information
    header = dlmread(filepath, ',', [0 0 0 0]);
    num_traj = num_traj + header(1,1);
    
    % Process the data
    % TODO - Is there a way to do this with only one file read?
    data = dlmread(filepath, ',', 1, 0);
    
    % Grab or validate the time vector
    if i == 1
        t = data(:,1);
        agg_data = zeros(size(data));
    else
        if norm(t - data(:,1)) ~= 0
            disp(['ERROR - Times in file ', files(i,1).name, ' do not align']);
            return;
        end
    end
    
    agg_data = agg_data + data;
end

% Normalize and unpack the aggregated data
agg_data = agg_data/num_traj;
[~, es, ess, es2, ess2] = unpack_symm(agg_data);
esSe = sqrt(es2 - es.^2)/sqrt(num_traj);
essSe = sqrt(real(ess2) - real(ess).^2)/sqrt(num_traj) + ...
     1i*sqrt(imag(ess2) - imag(ess).^2)/sqrt(num_traj);
 
 if min(min(es2 - es.^2)) < 0
     disp('WARNING: Found negative estimated variance in first order expectations');
 end
 
 if min(min(min(real(ess2) - real(ess).^2))) < 0 || ...
    min(min(min(imag(ess2) - imag(ess).^2))) < 0
    disp('WARNING: Found negative estimated variance in second order expectations');
 end