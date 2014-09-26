function [responseFns, regressor_events, tResponse, chanlocs] = LoadMultiLevelResponseFns(subjects, glmTypes)

% [responseFns, regressor_events, tResponse, chanlocs] =
%       LoadMultiLevelResponseFns(subjects, glmTypes)
% 
% INPUTS:
% - subjects is an N-element vector indicating the subject numbers.
% - glmTypes is an M-element cell array of strings indicating the suffixes
% of the GLM results you want to load. Ideally they should be a sequence
% (e.g., {'Saccade','Saccade-Type','Saccade-Type-SqNum'}).
%
% OUTPUTS:
% - responseFns is an M-element cell array of 4D matrices. responseFns{j} 
% is a set of response functions from all glms of glmType{j}. It is a
% (DxTxP(j)xN) matrix, where D=#channels, T=#timepoints, P(j)=#events in
% analysis glmType{j}.
% - regressor_events is an M-element cell array of cell arrays of strings.
% regressor_events{j} is a P(j)-element array containing the names of the 
% events in glmType{j}.
% - tResponse is a T-element vector indicating the times of each sample 
% relative to the locking event.
% - chanlocs is a D-element array of channel location structs taken from an
% eeglab file, which can be used for scalp maps using topoplot().
%
% Created 2/7/13 by DJ.
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 3/19/14 by DJ - tResponse in cells

nSubjects = numel(subjects);
nTypes = numel(glmTypes);
% glmTypes = {'LREvents-v2pt1','LREvents-NewType-v2pt1','LREvents-NewType-SqNum-v2pt1'};

%% Load Data
fprintf('Getting data from %d subjects...\n',nSubjects);
responseFns = cell(1,nTypes);
regressor_events = cell(1,nTypes);
for i=1:nSubjects
    fprintf('%d...',subjects(i))
    for j=1:nTypes
        S = load(sprintf('sq-%d-GLMresults-%s',subjects(i),glmTypes{j}));
        if isnumeric(S.responseFns)
            responseFns{j}(:,:,:,i) = S.responseFns;
        else
            responseFns{j}(:,:,:,i) = S.responseFns{S.iLevel};
        end
        regressor_events{j} = S.regressor_events{S.iLevel};
    end    
end
if isnumeric(S.tResponse)
    tResponse = S.tResponse;
else
    tResponse = S.tResponse{S.iLevel};
end
chanlocs = S.EEG.chanlocs;
disp('Done!')