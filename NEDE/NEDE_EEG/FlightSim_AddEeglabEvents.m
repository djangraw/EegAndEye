function [EEG, times, codes, iRemoved] = FlightSim_AddEeglabEvents(y,EEG,SYNC_CODES,eyeTimes,eyeCodes)

% Removes the current events from the EEGLAB file and replaces them with
% new ones according to the events_rule specified in the code.  Then
% re-saves the data, overwriting the old file.
%
% [EEG, times, codes, iRemoved] = FlightSim_AddEeglabEvents(x,EEG)
%
% NOTES:
% - EEGLAB should already be started.  
%
% INPUTS:
% -x
% -EEG
% -eyeTimes is a cell array of the times of events in each session (on the
% eyelink clock).
% -eyeCodes is a cell array of the corresponding event codes.
%
% OUTPUT:
% -EEG is the output file's eeglab struct.
%
% Created 9/25/14 by DJ based on NEDE_AddEeglabEvents.

%% CHECK INPUTS AND SET UP

tic

if ~exist('SYNC_CODES','var')
    SYNC_CODES = [];
end

% Get sync events
[eyeSyncs, eegSyncs] = NEDE_GetSyncEvents(y,EEG,SYNC_CODES);
if ischar(EEG.event(1).type)
    boundaryTimes = [0, EEG.event(strmatch('boundary', {EEG.event(:).type})).latency, EEG.pnts]*1000/EEG.srate; % in seconds
    % TEMPORARY KLUGE FIX FOR ONE SPECIFIC FILE
    if strncmp(EEG.filename,'FlightSim_1pt0-1-all',length('FlightSim_1pt0-1-all'))
        types =  {EEG.event(:).type};
        iSyncEvents = find(strcmp('208',types));
        times = [EEG.event(iSyncEvents).latency];
        iDelay = find(diff(times)>2000);
        boundaryTimes = unique([boundaryTimes, times(iDelay(3:end))]);    
    end
else
    boundaryTimes = [EEG.xmin, EEG.xmax];
end

if 0%~isequal(eegSyncs(:,2),eyeSyncs(:,2)) || ~isequal(eegSyncs(:,3),eyeSyncs(:,3)) % number of events that aren't in exactly the right spot
    error('Sync events don''t match up!');
end

nSessions = length(y);
[eegTimes,iRemoved] = deal(cell(1,nSessions));

%% Get events times and codes
for i=1:nSessions

    isInSession_sync = eyeSyncs(:,3)==i;
    eegTimes{i} = interp1(eyeSyncs(isInSession_sync,1),eegSyncs(isInSession_sync,1),double(eyeTimes{i}),'linear','extrap');
    iRemoved{i} = find(eegTimes{i}<boundaryTimes(i) | eegTimes{i}/EEG.srate>boundaryTimes(i+1));
    
    fprintf('Session %d: removing %d events\n',i,numel(iRemoved{i}));
    eegTimes{i}(iRemoved{i}) = [];
    if ~isempty(eyeCodes)
        eyeCodes{i}(iRemoved{i}) = [];
    end
end



%% Concatenate results
nEvents = numel(cat(1,eegTimes{:}));
times = reshape(cat(1,eegTimes{:}),nEvents,1)/1000; % get time in seconds
if ~isempty(eyeCodes)
    codes = reshape(cat(1,eyeCodes{:}),nEvents,1);
else
    codes = repmat({'event'},size(times));
end

%% sort results
[times,order] = sort(times,'ascend');
codes = codes(order);

%% Get events matrix and import into EEGLAB struct
if ~isempty(eyeCodes) % only actually import events if codes are specified.
    events = [num2cell(times), codes]; % the times (in s) and codes of each event
    assignin('base','events',events); % eeglab's importevent function grabs variable from base workspace
    EEG = pop_importevent( EEG, 'append','yes','event','events','fields',{'latency' 'type'},'timeunit',1,'optimalign','off');
    EEG = eeg_checkset( EEG );
end

toc