function visible_times = GetVisibleTimes(x,timekeeper)

% visible_times = GetVisibleTimes(x,timekeeper)
%
% - This function finds the times when each object entered or exited the 
% subject's view (a metric that we may use to anchor our EEG analyses).
% - INPUT x is the result of ImportData, which is usually stored in a .mat
% file called '3DS-[subject]-[session]'.
% - optional INPUT timekeeper should be 'eeg' or 'eyelink', and times will be
% reported in the units of eeg or eyelink samples (respectively). Default
% is 'eeg'.
% - OUTPUT visible_times is an nx3 matrix, where n is the number of times 
% an object was visible in the scene.  Each row is [object number, appear
% time, disappear time].
% Times are in units of EEG samples - divide by the EEG sampling rate to 
% get the time in seconds.
%
% Created 8/31/10 by DJ.
% Updated 8/4/11 by DJ - added date input to GetNumbers.

% Handle inputs
if nargin<2
    timekeeper = 'eeg';
end
if ~ismember({'eeg','eyelink'},timekeeper)
    error('Input ''timekeeper'' must be ''eeg'' or ''eyelink''!');
end

% SETUP
GetNumbers(x.date); % Get numbers struct
nObjects = numel(x.objects);
% Parse out stuff from data struct
event_times = x.(timekeeper).object_events(:,1);
event_codes = x.(timekeeper).object_events(:,2);
% Initialize output vectors
visible_times = [];

% FIND APPEARANCES OF EACH OBJECT
for i=1:nObjects
    % Find the relevant info
    appearTimes = event_times(event_codes==Numbers.ENTERS+i);
    disappearTimes = event_times(event_codes==Numbers.EXITS+i);
    if numel(appearTimes)<numel(disappearTimes)
        appearTimes = [appearTimes; x.eeg.start_time];
    elseif numel(disappearTimes)<numel(appearTimes)
        disappearTimes = [disappearTimes; x.eeg.end_time];
    end
    visible_times = [visible_times; repmat(i,numel(appearTimes),1), appearTimes, disappearTimes];
end

% SORT BY APPEAR TIME
[~,order] = sort(visible_times(:,2));
visible_times = visible_times(order,:);