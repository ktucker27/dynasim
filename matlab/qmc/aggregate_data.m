function [t, agg_data, num_traj] = aggregate_data(files)

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
agg_data(:,1) = t;