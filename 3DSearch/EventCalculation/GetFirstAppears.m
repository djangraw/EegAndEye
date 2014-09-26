function [targetApp, distractorApp] = GetFirstAppears(x,sec_between_repeats)

% [targetApp, distractorApp] = GetFirstAppears(x,sec_between_repeats)
%
% - This function finds the first time each object appeared in the 
% subject's view (a metric that we may use to anchor our EEG analyses).
% - INPUT x is the result of ImportData, which is usually stored in a .mat
% file called '3DS-[subject]-[session]'.
% - OUTPUTS targetApp and distractorApp are vectors of the first
% appearances of each target object and distractor object, respectively.
% They are in units of EEG samples - divide by the EEG sampling rate to get
% the time in seconds.
%
% Created 7/27/10 by DJ.
% Updated 7/29/10 by DJ (changed events field back to eyelink)
% Updated 11/5/10 by DJ (corrected line: appearTimes = [x.eeg.start_time...

% Handle inputs
if nargin==1
    sec_between_repeats = Inf;
end

% SETUP
GetNumbers; % Get numbers struct
nObjects = numel(x.objects);
% Convert sec_between_repeats to samples
fs = x.eeg.eventsamplerate;
time_between_repeats = sec_between_repeats*fs;
% Parse out stuff from data struct
event_times = x.eeg.object_events(:,1);
event_codes = x.eeg.object_events(:,2);
% Initialize output vectors
targetApp = [];
distractorApp = [];

% FIND FIRST APPEARANCE OF EACH OBJECT
for i=1:nObjects
    % Find the relevant info
    isTarget = strcmp(x.objects(i).tag,'TargetObject'); % Tag shouldbe 'TargetObject' or 'DistractorObject'
    appearTimes = event_times(event_codes==Numbers.ENTERS+i);
    disappearTimes = event_times(event_codes==Numbers.EXITS+i);
    if numel(appearTimes)<numel(disappearTimes)
        appearTimes = [x.eeg.start_time; appearTimes];
    elseif numel(disappearTimes)<numel(appearTimes)
        disappearTimes = [disappearTimes; x.eeg.end_time];
    end
    lastDisTime = -Inf; % last disappear time
    for j=1:numel(appearTimes)
        invisibleTime = appearTimes(j)-lastDisTime; % how long since we last saw this object?
        if invisibleTime < time_between_repeats % if it hasn't been long enough between repeats...
            lastDisTime = disappearTimes(j); % we treat this as if it were a continuation of the last appearance.
            continue; % skip this appearance
        else
            % Append to the appropriate vector
            if isTarget
                targetApp = [targetApp appearTimes(j)];
            else
                distractorApp = [distractorApp appearTimes(j)];
            end
            lastDisTime = disappearTimes(j);
        end
    end
end

% SORT BY TIME
targetApp = sort(targetApp);
distractorApp = sort(distractorApp);