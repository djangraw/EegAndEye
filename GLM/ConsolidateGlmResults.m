function ConsolidateGlmResults(prefix,subjects,suffixes,new_suffix)

% Combines multi-level analyses into one file to save disk space.
%
% ConsolidateGlmResults(prefix,subjects,suffixes,new_suffix)
%
% INPUTS:
% -prefix is a string indicating the experiment type
% -subjects is an N-element vector indicating the subject numbers
% -suffixes is an M-element cell array of strings indicating the filename
% suffixes for each subject's results (format is
% <prefix>-<subject>-GLMresults-<suffix>).
% -new_suffix is a string indicating the suffix filename you'd like to save
% the output as (same format as above).
%
% Created 4/30/13 by DJ.

% Declare constants
N = numel(subjects);
M = numel(suffixes);

% Main loop
for i=1:N    
    % Set up
    fprintf('---Subject %d ---\n',subjects(i));
    RF = cell(M,1);
    FN = cell(M,1);
    for j=1:M
        filename = sprintf('%s-%d-GLMresults-%s.mat',prefix,subjects(i),suffixes{j});
        % Load data
        fprintf('Loading %s...\n',filename)
        results = load(filename);
        % Extract important info
        if j==1
            EEG0 = results.EEG;
        end
        RF{j} = results.responseFns;
        FN{j} = filename;
    end
    % Overwrite two fields with 'consolidated' fields
    results.responseFns = RF;
    results.filenames = FN;
    results.EEG = EEG0;
    % Navigate into last file's folder
    fullpath = which(sprintf('%s-%d-GLMresults-%s.mat',prefix,subjects(i),suffixes{j}));
    folder = fileparts(fullpath);
    cd(folder)   
    % Save results
    fprintf('Saving %s-%d-GLMresults-%s.mat...\n',prefix,subjects(i),new_suffix);
    save(sprintf('%s-%d-GLMresults-%s.mat',prefix,subjects(i),new_suffix),'-struct','results');
    disp('Done!')
end
