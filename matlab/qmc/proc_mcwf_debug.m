function [t, jumps, state, num_traj] = proc_mcwf_debug(rootdir)

jumps = 0;
files = dir([rootdir, '/jumps*.txt']);
if size(files) > 0
    [t, jumps, num_traj] = aggregate_data(files);
end

state = 0;
files = dir([rootdir, '/state*.txt']);
if size(files) > 0
    [t2, state, num_traj2] = aggregate_data(files);
    
    if size(jumps) > 1
        if max(abs(t - t2)) ~= 0
            disp('WARNING: Times in debug directory do not match');
        end
        
        if num_traj ~= num_traj2
            disp('WARNING: Number of trajectories in debug directory do not match');
        end
    end
end
