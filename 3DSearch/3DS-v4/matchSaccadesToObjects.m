function ts_saccades = matchSaccadesToObjects(timestamps,ts_events)

% ts_saccades = matchSaccadesToObjects(timestamps,ts_events)
% -This function takes the results of find_writeioport and 
% find_saccadeevents and uses them to find which object the subject was
% looking at with each recorded saccade.
% -The INPUT timestamps is the result of find_saccadeevents.  ts_events is
% the result of find_writeioport.
% -The OUTPUT ts_saccades is an nx2 matrix, where n is the number of
% saccades in the experiment.  The first column is the eyelink timestamp of
% the end time of a saccade, and the second column is the number of the
% object to which a saccade was made.
% 
% Created 6/14/10 by DJ.
% Updated 6/17/10 by DJ - Made error checks
% Updated 7/14/10 by DJ - check for no timestamps

%Setup
GetNumbers;
ts_events = ts_events((ts_events(:,2)>Numbers.SACCADE_TO & ...
    ts_events(:,2)<Numbers.ENTERS),:); % crop to only the saccade events

if isempty(timestamps) % if there are no timestamps
    ts_saccades = []; % make sure output is defined
end

% Find the corresponding event for each timestamp
for i = 1:numel(timestamps)
    ts_saccades(i,1) = timestamps(i); % first column = timestamp
    event_index = find(ts_events(:,1)<=timestamps(i),1,'last'); % find most recent saccade events
    % check for errors
    if numel(event_index) ~=1
        warning('Event not found!');
        ts_saccades(i,2) = NaN;
    else
        ts_saccades(i,2) = ts_events(event_index,2); % second column = number of object saccaded to most recently
    end
end
