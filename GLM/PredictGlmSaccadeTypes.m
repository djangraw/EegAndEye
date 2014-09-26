function PredictGlmSaccadeTypes(EEG, responseFns, event_types, artifact_types, extent_ms, artifact_extent_ms)

% PredictGlmSaccadeTypes(EEG, responseFns, event_types, artifact_types, extent_ms, artifact_extent_ms)
%
% Created 1/30/12 by DJ (unfinished).
% Updated 3/22/13 by DJ - added extent inputs, use GetGlmRegressors_v2p0


disp('Setting up...');
% Set options

% Find constants
D = EEG.nbchan;
dt = 1000/EEG.srate;
regressor_range = round(extent_ms/dt); % how many samples should each response function extend?
artifact_range = round(artifact_extent_ms/dt); % how many samples should each artifact affect?
t = (1:EEG.pnts)*dt; % time vector for EEG

% Find event times
Nr = numel(event_types);
event_times = cell(1,Nr);
for i=1:Nr
    event_times{i} = [EEG.event(strcmp(event_types{i},{EEG.event(:).type})).latency]*dt;
end

% Find blink events
artifact_times = [EEG.event(ismember({EEG.event(:).type}, artifact_types)).latency]*dt;

disp('Getting regressors...');
% Get regressors
% [s,S] = GetGlmRegressors(t,event_times,artifact_times,Nt);
[s,S] = GetGlmRegressors_v2p0(t,event_times,artifact_times,regressor_range,artifact_range);

% Crop data
isRelevantTime = sum(S,2)>0;
trial_starts = find(diff(isRelevantTime)==1);
trial_ends = find(diff(isRelevantTime)==-1);
nTrials = length(trial_starts);

saccade_times = cell(1,nTrials);
for i=1:nTrials
    trial_time = trial_starts(i):trial_ends(i);
    saccade_times{i} = find(sum(s(:,trial_time)));
%     cla; hold on;
%     PlotVerticalLines(trial_time(saccade_times{i}),'r');
%     pause;
    Ns = numel(saccade_times{i});
    H = zeros(Nr,D*Nh*Ns);
    for j=1:numel(trial_time)
        trial_time(j);
    end
end
    