% GetGroupSvdResults_script.m
%
% A script to load GLM results from subjects 6-10, put them in the right
% format, and send them into GetGroupSvdResults.
%
% Created 4/17/12 by DJ.  It's specific for given subjects now, but can be modified.
% Updated 6/11/12 by DJ - added S14, weights/timcourses outputs

%% Set parameters
subjects = [9:11 13:15 17:19];
nComponents = 3;
nTrialTypes = 5;
timeWindow = [0 500];
nSubjects = numel(subjects);
glmType = 'LREvents-SqNum-NewType-v2pt1';


%% Load data (if you're sure the structs all match)
clear S
fprintf('Getting data from %d subjects...\n',nSubjects);
for i=1:nSubjects
    fprintf('%d...',subjects(i))
    S(i) =  load(sprintf('sq-%d-GLMresults-%s',subjects(i),glmType));    
end
disp('Done!')


%% Load Data
all = cell(1,nSubjects);
fprintf('Getting data from %d subjects...\n',nSubjects);
for i=1:nSubjects
    fprintf('%d...',subjects(i))
    if subjects(i) < 9
        all{i} = load(sprintf('sq-%d-GLMresults-Saccade-9reg',subjects(i)));                
    else
        all{i} = load(sprintf('sq-%d-GLMresults-%s',subjects(i),glmType));
    end    
end
disp('Done!')

%% Organize data
clear S
for i=1:nSubjects
%     if ~isfield(all{i},'nuisance_events')
%         all{i}.nuisance_events = {'Saccade'};        
%     end
    if ~isfield(all{i},'regressor_events')
        all{i}.regressor_events = all{i}.nuisance_events;
%         all{i}.nuisance_events = {};
    end
    S(i) = orderfields(all{i},all{1});
end
S = UpdateGlmResultsFormat(S);

%% Plot results
[weights,timecourse,avgWeights, avgTimecourses, chanlocs] = GetGroupSvdResults(S,nComponents,nTrialTypes,timeWindow);