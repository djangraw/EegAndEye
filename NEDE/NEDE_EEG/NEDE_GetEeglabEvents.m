function [eyeTimes, eyeCodes] = NEDE_GetEeglabEvents(y,EEG,SYNC_CODES)

% [eyeTimes, eyeCodes] = NEDE_GetEeglabEvents(y,EEG,SYNC_CODES)
%
% INPUTS:
% - y is an n-element NEDE data struct array.
% - EEG is an eeglab datastruct.
% - SYNC_CODES is a vector.
%
% OUTPUTS:
% - eyeTimes is an n-element cell vector in which eyeTimes{i} contains a
% p(i)-element vector of times when events occurred in session i.
% - eyeCodes is an n-element cell vector in which eyeCodes{i} contains a
% p(i)-element vector of the names of the events in session i.
%
% Created 9/26/14 by DJ based on NEDE_AddEeglabEvents.m.
% Updated 10/22/14 by DJ - added SYNC_CODES input, comments.

% Declare defaults
if nargin<3
    SYNC_CODES = []; % use default in NEDE_GetSyncEvents
end

% Get sync events
[eyeSyncs, eegSyncs] = NEDE_GetSyncEvents(y,EEG,SYNC_CODES);
if ischar(EEG.event(1).type)
    boundaryTimes = [0, EEG.event(strmatch('boundary', {EEG.event(:).type})).latency, EEG.pnts]*1000/EEG.srate; % in seconds
    % TEMPORARY KLUGE FIX FOR ONE SPECIFIC FILE
    if strncmp(EEG.filename,'FlightSim_1pt0-1-all',length('FlightSim_1pt0-1-all'))
        times = [EEG.event.latency];
        iDelay = find(diff(times)>2000);
        boundaryTimes = unique([boundaryTimes, times(iDelay(3:end))]);    
    end
else
    boundaryTimes = [EEG.xmin, EEG.xmax];
end

if ~isequal(eegSyncs(:,2),eyeSyncs(:,2)) || ~isequal(eegSyncs(:,3),eyeSyncs(:,3)) % number of events that aren't in exactly the right spot
    % TEMPORARY KLUGE FIX FOR ONE SPECIFIC FILE
    if strncmp(EEG.filename,'FlightSim_1pt0-1-all',length('FlightSim_1pt0-1-all'))
        warning('Sync events don''t match up... exception made for this file.');
    else
        error('Sync events don''t match up!');
    end
end


nSessions = length(y);
[eyeTimes,eegTimes,eyeCodes,iRemoved] = deal(cell(1,nSessions));

%% Get events times and codes
for i=1:nSessions
    x = y(i);
    % Parse inputs
    isTargetObject = strcmp('TargetObject',{x.objects(:).tag});
    % visibility events
    visible = x.events.visible;
    visible_objects = unique(visible.object);
    visible_isTarget = isTargetObject(visible_objects);
    visible_times = nan(length(visible_objects),2); 
    for j=1:numel(visible_objects)
        visible_times(j,:) = visible.time([...
            find(visible.object==visible_objects(j),1,'first'), ...
            find(visible.object==visible_objects(j),1,'last')]);
    end
    targetApp_times = visible_times(visible_isTarget,1);
    targetDisapp_times = visible_times(visible_isTarget,2);
    distApp_times = visible_times(~visible_isTarget,1);
    distDisapp_times = visible_times(~visible_isTarget,2);
    % saccade events
    saccade = x.events.saccade;
    saccadeToObj_times = saccade.time_end(saccade.isFirstToObject);
    saccadeToObj_obj = saccade.object_seen(saccade.isFirstToObject);
    saccadeToObj_isTarg = isTargetObject(saccadeToObj_obj);
    saccadeToTarg_times = saccadeToObj_times(saccadeToObj_isTarg);
    saccadeToDist_times = saccadeToObj_times(~saccadeToObj_isTarg);

    % construct matrix
    eyeTimes{i} = [targetApp_times; distApp_times; targetDisapp_times; distDisapp_times; ...
        saccadeToTarg_times; saccadeToDist_times; saccade.time_start; saccade.time_end; ...
        x.events.blink.time_start; x.events.blink.time_end; x.events.button.time; ...
        x.events.trial.time_start; x.events.trial.time_end];
    
    % Find event codes
    eyeCodes{i} = [repmat({'Targ Appear'},numel(targetApp_times),1); ...
        repmat({'Dist Appear'},numel(distApp_times),1); ...
        repmat({'Targ Disapp'},numel(targetDisapp_times),1); ...
        repmat({'Dist Disapp'},numel(distDisapp_times),1); ...
        repmat({'Targ Saccade'},numel(saccadeToTarg_times),1); ...
        repmat({'Dist Saccade'},numel(saccadeToDist_times),1); ...
        repmat({'Saccade Start'},numel(saccade.time_start),1); ...
        repmat({'Saccade End'},numel(saccade.time_end),1); ...
        repmat({'Blink Start'},numel(x.events.blink.time_start),1); ...
        repmat({'Blink End'},numel(x.events.blink.time_end),1); ...
        repmat({'Button Press'},numel(x.events.button.time),1); ...
        repmat({'Trial Start'},numel(x.events.trial.time_start),1); ...
        repmat({'Trial End'},numel(x.events.trial.time_end),1); ];

end
