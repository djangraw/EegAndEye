function EEG = AddBlinkMarkers(EEG0, y)

% EEG = AddBlinkMarkers(EEG0, y)
%
% INPUTS:
% -EEG0 is an eeglab struct made from a combined n sessions.
% -y is an n-element vector of behavior structs for the corresponding n 
% sessions, each including an 'eyelink' field with blink and saccade info.
%
% OUTPUTS:
% -EEG is the same as EEG0, but with 'BS' & 'BE' events for blink start & 
% end times.
%
% Created 6/4/13 by DJ.

% % Load data
% fprintf('loading data...\n');
% eegFilename = sprintf('%s-%d-all%s.set',prefix,subject,eegSuffix);
% EEG = pop_loadset(eegFilename);
% y = loadBehaviorData(subject,sessions,prefix);

% Declare variables
nSessions = numel(y);
allEvents = cell(1,nSessions);
allCodes = cell(1,nSessions);

% Get HEOG saccades
for i=1:nSessions
    if isfield(y(i),'eyelink')
        % Get blink-containing saccades
        iBlink = FindBlinkSaccades(y(i).eyelink.blink_times, y(i).eyelink.saccade_start_times, y(i).eyelink.saccade_times);

        % Get saccade times       
        startTimes = y(i).eyelink.saccade_start_times(iBlink);
        endTimes = y(i).eyelink.saccade_times(iBlink);
    else
        % Get blink-containing saccades
        iBlink = FindBlinkSaccades(y(i).blink.start_time, y(i).saccade.start_time, y(i).saccade.end_time);

        % Get saccade times       
        startTimes = y(i).saccade.start_time(iBlink);
        endTimes = y(i).saccade.end_time(iBlink);
    end
    % Make events matrix
    allEvents{i} = cat(1,startTimes,endTimes);
    allCodes{i} = cat(1,repmat({'BS'},size(startTimes)),...
        repmat({'BE'},size(endTimes)));
    
end

% Add events to EEG struct
y = AddSyncField(y);
EEG = AddEeglabEvents_MultiSession(EEG0,y,allEvents,allCodes);