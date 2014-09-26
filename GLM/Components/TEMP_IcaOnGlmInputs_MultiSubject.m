% TEMP_IcaOnGlmInputs_MultiSubject.m
% originally called TEMP_GetGlmIca_MultiSubject.m
%
% Created 2/1/13 by DJ for one-time use.
% Updated 4/24/13 by DJ - added load & plot, comments.


%% Set up
subjects = [9:11, 13:15, 17:19];
folders = {'2012-03-29-Pilot', '2012-04-09-Pilot', '2012-04-18-Pilot', ...
    '2012-05-04-Pilot', '2012-05-29-Pilot', '2012-06-12-Pilot', ...
    '2012-06-15-Pilot', '2012-06-21-Pilot', '2012-06-22-Pilot'};
%% Run
for i=1:numel(subjects)
    fprintf('--- SUBJECT %d ---\n',subjects(i))
    cd(folders{i});
    % Apply ICA to the GLM inputs
    [icaweights, icasphere, icawinv] = ApplyIcaToGlmInput(R(i));
    % Save results
    save(sprintf('sq-%d-GlmIca-LREvents-SqNum-NewType-v2pt1',subjects(i)),'icaweights','icasphere','icawinv');
    cd ..
end



%% Load & Plot
% Use TEMP_IcaOnGlmInputs_PlotResults.
