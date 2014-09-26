function [h, y, iEvents, truth] = PrepareGlmDataForClassification(results,EEG,iLevel_input)

% Created 1/31/12 by DJ.
% Updated 2/1/12 by DJ - changed name (was PreparedGlmDataForBryan)
% Updated 3/22/13 by DJ - use asymmetric influence, GetGlmRegressors_v2p0
% Updated 4/30/13 by DJ - responseFns in cells

if nargin<3
    iLevel_input = [];
end

disp('Setting up...');
% Set options
results = UpdateGlmResultsFormat(results); % new 3/22/13
UnpackStruct(results);
if ~isempty(iLevel_input)
    iLevel = iLevel_input;
end

% Find constants
dt = 1000/EEG.srate;
regressor_range = round(influence/dt); % how many samples should each response function extend?
artifact_range = round(artifact_influence/dt); % how many samples should each artifact affect?
t = (1:EEG.pnts)*dt; % time vector for EEG

% Find event times
Nr = numel(regressor_events{iLevel});
event_times = cell(1,Nr);
for i=1:Nr
    event_times{i} = [EEG.event(strcmp(regressor_events{iLevel}{i},{EEG.event(:).type})).latency]*dt;
end

% Find blink events
artifact_times = [EEG.event(ismember({EEG.event(:).type}, artifact_events)).latency]*dt;

disp('Getting regressors...');
% Get regressors
% [s,S] = GetGlmRegressors(t,event_times,artifact_times,Nt);
[s,S] = GetGlmRegressors_v2pt0(t,event_times,artifact_times,regressor_range,artifact_range);

disp('Preparing output...');
% crop data
isRelevantTime = sum(S,2)>0;
sCrop = s(:,isRelevantTime);
y = EEG.data(:,isRelevantTime);

[truth, iTimeZero] = find(sCrop);
iEvents = iTimeZero + regressor_range(1);
h = permute(responseFns{iLevel},[3 2 1]);