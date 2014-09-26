% ExpandLastGlmLevel_script
%
% Created 3/19/14 by DJ.

prefix = 'sf';
oldSuffix = 'Events-SqNum-Type-v3pt1.mat';
newSuffix = 'Events-SqNum-Type-v3pt2.mat';

% --- LOAD OLD DATA
[R,subjects,basedir,folders] = LoadAllGlmResults(prefix,oldSuffix);
R = UpdateGlmResultsFormat(R);
%%
% --- REDO DESIRED LEVELS OF GLM
cd(basedir)

for i=1:numel(subjects)
    % SET UP
    fprintf('--- SUBJECT %d ---\n',subjects(i));
    cd(folders{i});
    subject = subjects(i);
    
    % UPDATE PARAMETERS
    clear params
    params.filenames = R(i).filenames;
    params.filenames{end} = sprintf('%s-%d-GLMresults-%s',prefix,subject,newSuffix);
    params.influence = R(i).influence;
    params.influence{end} = [0 750];
    
    % REDO ANALYSIS
    if ~exist(params.filenames{end},'file')
        RedoGlmLevels(R(i),params);
    else
        disp('--- SKIPPING REGULAR GLM ---');
    end
    
    cd ..
    fprintf('--- DONE! ---\n');
end
