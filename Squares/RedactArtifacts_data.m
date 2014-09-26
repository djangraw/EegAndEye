function EEG = RedactArtifacts_data(EEG,artifact_events,tRange)

% Replace the data around artifacts with NaN's.
%
% RedactArtifacts_data(EEG,artifact_events,tRange)
%
% INPUTS:
% -EEG is an eeglab data struct (not epoched) with events already added.
% -artifact_events is a string or cell array of strings indicating the
% events.
% -tRange is a 2-element vector indicating the minimum and maximum time (in
% ms, relative to the event time) that should be redacted.
%
% Created 3/20/13 by DJ.

if ~iscell(artifact_events)
    artifact_events = {artifact_events};
end
if size(EEG.data,3)>1
    error('This function is for use with continuous datasets only!');
end


sampleRange = [ceil(tRange(1)/1000*EEG.srate), floor(tRange(2)/1000*EEG.srate)]; 

iEvents = find(ismember({EEG.event.type},artifact_events));

for i=1:numel(iEvents) % for each event
    sampleEvent = round(EEG.event(iEvents(i)).latency); % time when event occurred
    samplesToRedact = max(sampleEvent+sampleRange(1), 1):min(sampleEvent+sampleRange(2), EEG.pnts); % samples to redact (must fall within limits of data)
    EEG.data(:,samplesToRedact) = NaN; % REDACTED!
end