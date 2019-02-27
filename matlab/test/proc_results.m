function reports = proc_results(rootdir, do_plots)

% proc_results: Processes a directory tree of output from multiple
%               simulators comparing results
%
% INPUT: rootdir  - Root of the results directory structure
%        do_plots - Will produce plots of spin squeezing for all results
%                   if provided and non-zero
%
% OUTPUT: TODO

PLOT_ROWS = 2;
PLOT_COLS = 3;
PLOTS_PER_PAGE = PLOT_ROWS*PLOT_COLS;

if nargin < 2
    do_plots = 0;
end

if ~strcmp(rootdir(1,size(rootdir,2)), '/')
    rootdir = [rootdir, '/'];
end

runTypes = {'symm'; 'mcwf'};
fun = @(name)(strcmp(name, '.'));

% Populate results directory map
resdirs = containers.Map;
reports = containers.Map;
for i = 1:size(runTypes,1)
    files = dir([rootdir, '/**/', runTypes{i}, '/']);
    files = files(cellfun(fun, {files.name}));
    for fileIdx = 1:size(files,1)
        % Strip off rootdir and type from results directory to use as key
        resdir = strrep(files(fileIdx).folder, rootdir, '');
        reslist = strsplit(resdir, '/');
        resdir = strjoin(reslist(1, 1:size(reslist, 2) - 1), '/');
        
        % Initialize structs if this is the first time we've seen this key
        if ~resdirs.isKey(resdir)
            resdirs(resdir) = [];
            
            reports(resdir) = new_report();
        end
        
        % Add run type to file info
        files(fileIdx).type = runTypes{i};
        
        % Add file info to results map
        resdirs(resdir) = [resdirs(resdir); files(fileIdx)];
        
        % Add this run type to the report
        report = reports(resdir);
        report.types = [report.types, runTypes{i}];
        reports(resdir) = report;
    end
end

% Iterate over result directories processing results
reskeys = resdirs.keys;
for resIdx = 1:size(reskeys, 2)
    resdir = reskeys{resIdx};
    dirlist = resdirs(resdir);
    report = reports(resdir);
    
    disp(['Processing directory: ', resdir]);
    
    types = {};
    for fileIdx=1:size(dirlist,1)
        fileInfo = dirlist(fileIdx);
        
        types = [types, fileInfo.type];
        
        disp(['Found type: ', fileInfo.type]);
        
        % Process the results directory
        if fileInfo.type == 'symm'
            results(fileIdx) = proc_symm(fileInfo.folder);
        elseif fileInfo.type == 'mcwf'
            results(fileIdx) = proc_mcwf(fileInfo.folder);
        else
            disp(['ERROR: Unknown file type ', fileInfo.type]);
            return;
        end
    
        if results(fileIdx).init == 0
            disp('Error processing results');
            return;
        end
        
        % Compare current results against all previous types
        for tidx = 1:fileIdx-1
            sq2db1 = 10*log10(real(results(tidx).sq2));
            sq2db2 = 10*log10(real(results(fileIdx).sq2));
            report.sq2_uniform(tidx, fileIdx-1) = interp_compare(results(tidx).t, ...
                                                                 sq2db1, ...
                                                                 results(fileIdx).t, ...
                                                                 sq2db2);
        end
        
        % Plot the results if requested
        if do_plots
            if fileIdx == 1
                if mod(resIdx-1, PLOTS_PER_PAGE) + 1 == 1
                    figure('Color',[1,1,1],'Position', [1,1,1000,800]);
                end
                subplot(PLOT_ROWS, PLOT_COLS, mod(resIdx-1, PLOTS_PER_PAGE) + 1);
                hold on
                grid on
                title(resdir);
            end
            
            plot(results(fileIdx).t, 10*log10(real(results(fileIdx).sq2)));
            lgd = legend(types);
            lgd.FontSize = 12;
        end
    end
    
    % Save the modified report
    % TODO - Is there a better way to do this?
    reports(resdir) = report;
end
end

function n = extract_n(filepath)
n = 0;
pat = '.*/N([^/]+)/.*';
tok = regexp(filepath, pat, 'tokens');
if size(tok, 1) ~= 1
    disp(['ERROR: Found incorrect number of N ', num2str(size(tok,1)), ' in filepath ', filepath]);
    return;
end
n = str2double(tok{1});
end

function results = proc_symm(resdir)
results.init = 0;
results.t = [];
results.sq2 = [];
files = dir([resdir, '/symm_N*.txt']);
if size(files, 1) ~= 1
    disp(['ERROR: Found incorrect number of symm files: ', num2str(size(files,1))]);
    return;
end

n = extract_n(files(1).folder);
if n == 0
    return;
end

C = dlmread([files(1).folder, '/', files(1).name]);
[results.t, results.es, results.ess] = unpack_symm(C);
results.sq2 = get_squeezing(n, results.es, results.ess);
results.init = 1;
end

function results = proc_mcwf(resdir)
results.init = 0;
results.t = [];
results.sq2 = [];

n = extract_n(resdir);
if n == 0
    return;
end

[results.t, results.es, results.ess, ~] = proc_mcwf_output(resdir, 0, 1);
results.sq2 = get_squeezing(n, results.es, results.ess);
results.init = 1;
end

function report = new_report()
report = struct();
report.types = {};
end

function idx = type_idx(runTypes, type)
idx = -1;
for i=1:size(runTypes,1)
    if strcmp(runTypes{i}, type)
        idx = i;
        break;
    end
end
end