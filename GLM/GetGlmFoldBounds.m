function foldBounds = GetGlmFoldBounds(EEG,regressor_events,artifact_events,extent_ms,nFolds,artifact_extent_ms)

% Separates a squares expt EEG file into N separate structs.
%
% foldBounds = GetGlmFoldBounds(EEG,regressor_events,artifact_events,
%  extent_ms,nFolds,artifact_extent_ms)
%
% INPUTS:
% - EEG is an eeglab data struct from a squares experiment.
% - regressor_events, artifact_events, extent_ms, and artifact_extent_ms 
% are used in a way that parallels the other GLM programs.
% - nFolds is a scalar indicating the number of structs you would like to
% create.
%
% OUTPUTS:
% - foldBounds is an nFolds-element vector in which the bounds of data used
% in fold i will be between samples foldBounds(i) and foldBounds(i+1).
%
% Is careful to recommend separations between trials.  Usually called by
% RunXfoldEegGlm().
%
% Created 2/1/12 by DJ.
% Updated 3/19/12 by DJ - comments
% Updated 3/22/13 by DJ - use asymmetric influence, GetGlmRegressors_v2p0

if nargin<6 || isempty(artifact_extent_ms)
    warning('artifact_extent_ms not specified!');
    artifact_extent_ms = extent_ms;
end

% Find constants
dt = 1000/EEG.srate;
% Nt = round(extent_ms/dt); % how many samples should each response function extend?
regressor_range = round(extent_ms/dt); % how many samples should each response function extend?
artifact_range = round(artifact_extent_ms/dt); % how many samples should each artifact affect?
t = (1:EEG.pnts)*dt; % time vector for EEG

% Find event times
Nr = numel(regressor_events);
event_times = cell(1,Nr);
for i=1:Nr
    event_times{i} = [EEG.event(strcmp(regressor_events{i},{EEG.event(:).type})).latency]*dt;
end

% Find blink events
artifact_times = [EEG.event(ismember({EEG.event(:).type}, artifact_events)).latency]*dt;

disp('Getting regressors...');
% Get regressors
% [s,S] = GetGlmRegressors(t,event_times,artifact_times,Nt);
[s,S] = GetGlmRegressors_v2p0(t,event_times,artifact_times,regressor_range,artifact_range);

% Get info about data
isRelevantTime = sum(S,2)>0; % Will this sample be used in the final GLM?
nEventsSoFar = cumsum(sum(s,1)); % How many regressor events have taken place so far?
nEventsPerFold = sum(s(:)/nFolds); % How many events should be included per fold?

% Translate this info into fold bounds
foldBounds = nan(1,nFolds+1);
foldBounds(1) = 1;
foldBounds(end) = EEG.pnts;
for i=1:nFolds-1    
    foldBounds(i+1) = find(nEventsSoFar>=nEventsPerFold*i & ~isRelevantTime',1);    
end
disp('Done!')