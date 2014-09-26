function [trialtimes trialtypes] = get_trialtimes_portversion(ts_events)

% trialtimes = get_trialtimes_portversion(ts_events)
% -Given the timestamps and events from the eyetracker, finds and returns 
% the time when each trial was loaded, started, and ended.  
% -Input ts_events is an nx2 matrix, where n is the number of events.  The 
% first column is the timestamp of each event, the second column is the
% numeric code of that event (see script GetNumbers to decode them).
% -Output trialtimes is an mx3 matrix, where m is the number of trials.
% Column 1 is the timestamp of the trial loading, column 2 is the timestamp
% of the trial starting, and column 3 is the timestamp of the trial ending.
% -This program is obsolete, but we're keeping it because it might be
% useful in a future version of 3DSearch.
%
% Created 6/15/10 by DJ.
% Updated 6/17/10 by DJ - added trialtypes output and error check
% Updated 7/26/10 by DJ - just changed the name.

GetNumbers; % get constants from currrent Unity experiment

% Get the times
load_times = ts_events(ts_events(:,2)==Numbers.LOAD_TRIAL,1); % all the load_trial times in a vector
start_times = ts_events(ismember(ts_events(:,2),...
    Numbers.START_TRIAL + [Numbers.STATIONARY, Numbers.MOVING, ... 
    Numbers.POPUP]),1); % all the start_trial times in a vector
trialtypes = ts_events(ismember(ts_events(:,2),...
    Numbers.START_TRIAL + [Numbers.STATIONARY, Numbers.MOVING, ... 
    Numbers.POPUP]),2) - Numbers.START_TRIAL; % the surpriseLevel variable for each trial
end_times = ts_events(ts_events(:,2)==Numbers.END_TRIAL,1); % all the end_trial times in a vector

% Make sure they're the same length
if ~isequal(numel(load_times),numel(start_times),numel(end_times))
    error('Numbers of load, start, and end trial codes don''t match!');
end

% Create the final product
trialtimes(:,1) = load_times;
trialtimes(:,2) = start_times;
trialtimes(:,3) = end_times;