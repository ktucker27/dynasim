function [t, es, ess, num_traj, esSe, essSe, jumps, state] = proc_mcwf_output(rootdir, proc_dbg, quiet)

% proc_mcwf_output: Aggregates and organizes results coming out of
%                   the MCWF Java executable
%
% INPUT: rootdir  - Can be a directory of files to aggregate, or a
%                   single file to process
%        proc_dbg - Process debug output if nonzero
%           quiet - Supress command line output if nonzero
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

if nargin < 3
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

[t, agg_data, num_traj] = aggregate_data(files, quiet);

%dlmwrite('/Users/tuckerkj/output/20181112/test/agg.csv', agg_data, 'delimiter', ',', 'precision', 12);

[~, es, ess, es2, ess2] = unpack_symm(agg_data);

% Calculate statistical error
esSe = sqrt(es2 - es.^2)/sqrt(num_traj);
essSe = sqrt(real(ess2) - real(ess).^2)/sqrt(num_traj) + ...
     1i*sqrt(imag(ess2) - imag(ess).^2)/sqrt(num_traj);
 
% Process debug directory if found
dbgdir = [rootdir, '/debug'];
jumps = 0;
state = 0;
if nargin > 1 && proc_dbg ~= 0 && isfolder(dbgdir)
    if ~quiet
        disp('Found debug directory, aggregating results...');
    end
    [t2, jumps, state, num_traj2] = proc_mcwf_debug(dbgdir, quiet);
    if max(abs(t - t2)) ~= 0
        disp('WARNING: Difference found between output and debug times');
    end
    
    if num_traj ~= num_traj2
        disp('WARNING: Difference found between output and debug number of trajectories');
    end
end

if min(min(es2 - es.^2)) < -1.0e-7
    disp(['WARNING: Found negative estimated variance in first order expectations: ', num2str(min(min(es2 - es.^2)))]);
end

if min(min(min(real(ess2) - real(ess).^2))) < -1.0e-7 || ...
        min(min(min(imag(ess2) - imag(ess).^2))) < -1.0e-7
    disp(['WARNING: Found negative estimated variance in second order expectations: ', num2str(min(min(min(real(ess2) - real(ess).^2)))), ' + i', num2str(min(min(min(imag(ess2) - imag(ess).^2))))]);
end