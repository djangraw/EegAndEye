function [EEG, times, codes, iRemoved] = AddEeglabEvents_MultiSession(EEG, y, elTimes, elCodes)

% [EEG, times, codes, iRemoved] = AddEeglabEvents_MultiSession(EEG, y, elTimes, elCodes)
%
% INPUTS:
% -EEG is an EEGLAB data struct containing multiple concatenated data files
% -y is an n-element vector of squares structs.
% -elTimes is an n-element column vector of cells in which cell i contains
% the eyelink times in session i.
% -elCodes is an n-element column vector of cells in which cell i contains 
% a cell array of event type strings for session i.
%
% OUTPUTS:
% -EEG is the input struct with the requested event times added.
% -times is a column vector of the EEG times of the specified events.
% -codes is a column vector of cells containing the string codes for the
% specified events.
% -iRemoved is a column vector of cells containing the indices of elTimes
% and elCodes that fell outside the legal area and so were removed.
%
% Created 11/28/11 by DJ.
% Updated 5/9/12 by DJ - added times/codes outputs, default for elCodes
% Updated 4/23/13 by DJ - added iRemoved output to avoid removing events
% with eeg_checkset

%% Set up
if nargin<4
    elCodes = [];
end
nSessions = numel(elTimes);
boundaryTimes = [0 EEG.event(strmatch('boundary', {EEG.event(:).type})).latency]/EEG.srate; % in seconds
if numel(boundaryTimes)<nSessions-1
    error('%d sessions given as input but only %d boundary events!',...
        nSessions,numel(boundaryTimes))
end

%% Convert to EEG times
eegTimes = cell(size(elTimes));
iRemoved = cell(size(elTimes));
for i=1:nSessions
    % Set up syncs
    elSyncs = y(i).sync.eyelink;
    eegSyncs = y(i).sync.eeg/y(i).eeg.eventsamplerate + boundaryTimes(i); % in seconds
    % Convert eyelink to eeg times
    eegTimes{i} = EyelinkToEegTimes(elTimes{i},elSyncs,eegSyncs);
    % Remove any events that fall outside the legal range
    iRemoved{i} = find(eegTimes{i}<EEG.xmin | eegTimes{i}>EEG.xmax);
    fprintf('Session %d: removing %d events\n',i,numel(iRemoved{i}));
    eegTimes{i}(iRemoved{i}) = [];
    if ~isempty(elCodes)
        elCodes{i}(iRemoved{i}) = [];
    end
end

%% Concatenate results
nEvents = numel(cat(1,eegTimes{:}));
times = reshape(cat(1,eegTimes{:}),nEvents,1); % get time in seconds
if ~isempty(elCodes)
    codes = reshape(cat(1,elCodes{:}),nEvents,1);
else
    codes = repmat({'event'},size(times));
end

%% Get events matrix and import into EEGLAB struct
if ~isempty(elCodes) % only actually import events if codes are specified.
    events = [num2cell(times), codes]; % the times (in s) and codes of each event
    assignin('base','events',events); % eeglab's importevent function grabs variable from base workspace
    EEG = pop_importevent( EEG, 'append','yes','event','events','fields',{'latency' 'type'},'timeunit',1,'optimalign','off');
    EEG = eeg_checkset( EEG );
end