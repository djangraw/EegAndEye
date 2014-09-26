% RedoWithRidgeRegression_script.m
%
% Created 8/14/14 by DJ.

cd(basedir)
prefix = 'sf';
subjects = [1:10 12:13];
old_suffix = 'Type-v3pt6-RampUp.mat';

new_params = struct('method','ridge');
new_params.lambda = {'auto'};
new_params.demeanX = false;
new_suffix = 'Type-v3pt6-RampUp-ridgeauto-nodemean.mat';
% new_params = struct('method','leastsquares','lambda',0);
% new_suffix = 'Type-v3pt6-RampUp.mat';

for i=1:numel(subjects)
    fprintf('--- SUBJECT %d ---\n',subjects(i));
    cd(folders{i});
    subject = subjects(i);

    if ~exist(sprintf('%s-%d-GLMresults-%s',prefix,subject,new_suffix),'file')
        R = load(sprintf('%s-%d-GLMresults-%s',prefix,subject,old_suffix));
        % NOTE: THIS WORKS FOR ONE-LEVEL GLM ONLY!        
        new_params.filenames = { sprintf('%s-%d-GLMresults-%s',prefix,subject,new_suffix) };       
        RedoGlmLevels(R,new_params);        
    else
        disp('--- SKIPPING REGULAR GLM ---');
    end
    cd ..
    fprintf('--- DONE! ---\n');
end

%% FIX SQUARES

% fake_filename = 'sf-13-GLMresults-Type-v3pt6-RampUp-ridge100.mat';
% new_suffix = 'Type-v3pt6-RampUp.mat';
% cd(basedir)
% 
% for i=1:numel(subjects)
%     fprintf('--- SUBJECT %d ---\n',subjects(i));
%     cd(folders{i});
%     subject = subjects(i);
% 
%     if exist([cd '/' fake_filename],'file')
%         delete([cd '/' fake_filename]);
%         load([cd '/' fake_filename]);
%         % NOTE: THIS WORKS FOR ONE-LEVEL GLM ONLY!        
%         filenames = { sprintf('%s-%d-GLMresults-%s',prefix,subject,new_suffix) };       
%         lambda = 0;
%         fprintf('Saving %s...',filenames{iLevel});
%         save(filenames{iLevel},'EEG','responseFns','tResponse','regressor_events',...
%             'filenames','iLevel','artifact_events','dataset','offset','influence',...
%             'artifact_influence','stddev','vthresh','method','trial_rej_rules','lambda');
%     else
%         disp('--- SKIPPING REGULAR GLM ---');
%     end
%     cd ..
%     fprintf('--- DONE! ---\n');
% end