function [targetSac, distractorSac] = GetFirstSaccades(x,sec_between_repeats)

% [targetSac, distractorSac] = GetFirstSaccades(x,sec_between_repeats)
%
% - This function finds the first saccade to each object (a metric that we
% may want to use to anchor our EEG analyses).
% - INPUT x is the result of ImportData, which is usually stored in a .mat
% file called '3DS-[subject]-[session]'.
% - INPUT sec_between_repeats is a the time (in s) an object must be out
% of view before we allow it to be considered a 'new' object (default: Inf)
% - OUTPUTS targetSac and distractorSac are vectors of the first saccades
% the subject made to each target and distractor object, respectively.
% They are in units of EEG samples - divide by the EEG sampling rate to get
% the time in seconds.
% 
% Created 7/27/10 by DJ.
% Updated 7/29/10 by DJ (updated, added sec_between_repeats input)
% Updated 8/19/10 by DJ (saccade_events and object_events are now separate)
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
saccade_times = x.eeg.saccade_events(:,1);
saccade_codes = x.eeg.saccade_events(:,2);
% Initialize output vectors
targetSac = [];
distractorSac = [];

% FIND FIRST SACCADE TO EACH OBJECT
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
        firstSac = saccade_times(find(saccade_codes==Numbers.SACCADE_TO+i & saccade_times>appearTimes(j) & saccade_times<disappearTimes(j),1)); % Find first saccade during this appearance
        if isempty(firstSac)
            continue; % do not even include this as an appearance - if we had a saccade in an appearance right after this one, we'd want to count it.
        elseif invisibleTime < time_between_repeats % if it hasn't been long enough between repeats, or if no saccade was made...
            lastDisTime = disappearTimes(j); % include this as an appearance - if we had a saccade in an appearance right after this one, we WOULDN'T want to count it.
            continue; % skip this appearance
        else
            % Append to the appropriate vector
            if isTarget
                targetSac = [targetSac firstSac];
            else
                distractorSac = [distractorSac firstSac];
            end
            % Update the time when this object was last seen and disappeared.
            lastDisTime = disappearTimes(j);
        end
    end
end

% SORT BY TIME
targetSac = sort(targetSac);
distractorSac = sort(distractorSac);