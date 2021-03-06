function EEG = RemoveEyeBlinkTrials(EEG,bad_event,t_window_epoch,t_window_saccade,t_discrim,experimentType)

% EEG = RemoveEyeBlinkTrials(EEG,t_window,bad_event)
%
% Remove any trials with eye blinks (or other 'bad events') in the
% specified time interval.
%
% INPUTS:
% - EEG is an EEGLAB data structure.
% - bad_event is a scalar specifying the number of the event to be removed
%   (see GetNumbers for details).
% - t_window_epoch is a 2-element vector specifying times (in ms) of the 
%   epoch in which blinks will affect a trial (and so should be removed). 
%   (default: [-200 0] to avoid baseline)
% - t_window_saccade is a 2-element vector specifying times (in ms) relative
%   to the first and last saccade in an epoch in which blinks will affect a 
%   trial (and so should be removed). (default: [-500 500] to stay outside
%   of TrainingWindowOffsets)
% - t_discrim is a vector containing the time range (in ms) relative to
%   t=0 in which saccades will be considered potential discriminating
%   saccades. (default: [0 1000], see computeSaccadeJitterPrior)
%
% OUTPUTS:
% - EEG is the inputted eeglab data structure, with eye blink trials removed
%   as specified.
%
% Created 12/1/10 by DJ
% Updated 12/7/10 by DJ - comments
% Updated 2/28/11 by DJ - made a function. 
% Updated 7/27/11 by DJ - added t_window_saccade input, changed defaults
% Updated 8/1/11 by DJ - consider both start and end saccade times, add
%    t_discrim input
% Updated 8/5/11 by DJ - updated defaults
% Updated 11/1/11 by DJ - added experimentType input to work with Squares
%    experiments

%% Set up
GetNumbers; % load Numbers Struct
Constants = GetSquaresConstants;
if nargin<2 || isempty(bad_event)    
    bad_event = Numbers.BLINK;
end
if nargin<3 
    t_window_epoch = [-200 0]; % time window (in ms) relative to t=0 in which to check for blinks
end
if nargin<4 
    t_window_saccade = [-500 500]; % time window (in ms) relative to saccades in which to check for blinks
end
if nargin<5
    t_discrim = [0 1000]; % epoch time range (min, max in ms) in which saccades could be used in discrimination.
end
if nargin<6 || isempty(experimentType)
    experimentType = '3DSearch';
end


%% Apply threshold and reject offending trials
fprintf('---Excluding trials with blinks at t=[%d, %d] or [%d, %d] relative to saccade\n',...
    t_window_epoch(1),t_window_epoch(2), t_window_saccade(1),t_window_saccade(2))
% saccadeTimes = GetEpochSaccades(EEG,{'SACCADE_START','SACCADE_END'},t_discrim);


%% Find eye blink trials
isBlinky = zeros(1,numel(EEG.epoch)); % was there a blink in the time window on this trial?
for i=1:numel(EEG.epoch)
    % 0. Extract necessary info from EEG struct
    % Get events in this epoch
    eventTypes = str2double([EEG.epoch(i).eventtype]);
    eventLatencies = cell2mat([EEG.epoch(i).eventlatency]);
    % Separate out events
    blinks = eventLatencies(eventTypes==bad_event);
    if strcmp(experimentType,'3DSearch')
        saccadeStart = eventLatencies(eventTypes==Numbers.SACCADE_START);
        saccadeEnd = eventLatencies(eventTypes==Numbers.SACCADE_END);
    elseif strcmp(experimentType,'Squares')
        saccadeStart = eventLatencies(ismember(eventTypes,Constants.SACCADESTART_BASE:Constants.SACCADEEND_BASE));
        saccadeEnd = eventLatencies(ismember(eventTypes,Constants.SACCADEEND_BASE:(2*Constants.SACCADEEND_BASE-Constants.SACCADESTART_BASE)));
    end
    % Get real start and end of each blink (defined by the saccade around it)
    blinkStart = blinks;
    blinkEnd = blinks;
    for j=1:numel(blinks)
        iBlinkStart = find(saccadeStart<blinks(j),1,'last');
        if ~isempty(iBlinkStart)
            blinkStart(j) = saccadeStart(iBlinkStart);
        end
        iBlinkEnd = find(saccadeEnd>blinks(j),1,'first');
        if ~isempty(iBlinkEnd)
            blinkEnd(j) = saccadeEnd(iBlinkEnd);
        end
    end
    
    
    % 1. Exclude trials with blinks in the specified epoch time window
    if  any(blinkStart>t_window_epoch(1) & blinkStart<t_window_epoch(2)) || any(blinkEnd>t_window_epoch(1) & blinkEnd<t_window_epoch(2))
        isBlinky(i) = 1; % Blink trial!
        continue;
    end
    
    
    % 2. Exclude trials with no saccades in the discrimination window
    % Find saccades within the discrimination window
    okSaccadeStart = saccadeStart(saccadeStart>=min(t_discrim) & saccadeStart<=max(t_discrim));
    okSaccadeEnd = saccadeEnd(saccadeEnd>=min(t_discrim) & saccadeEnd<=max(t_discrim)); 
    % For start and end times (separately), make sure you have at least one
    if isempty(okSaccadeStart)
        okSaccadeStart = saccadeStart(find(saccadeStart>=min(t_discrim),1,'first'));
        if isempty(okSaccadeStart)
            isBlinky(i)=1; % No saccade starts in this trial... exclude it.
            continue;
        end
    end
    if isempty(okSaccadeEnd)
        okSaccadeEnd = saccadeEnd(find(saccadeEnd>=min(t_discrim),1,'first'));
        if isempty(okSaccadeEnd)
            isBlinky(i)=1; % No saccade ends in this trial... exclude it.
            continue;
        end
    end
    
    % 3. Exclude trials with blinks in the specified saccade-locked time window
    % Get first and last saccades
    firstSaccade = min([okSaccadeStart, okSaccadeEnd]);  
    lastSaccade = max([okSaccadeStart, okSaccadeEnd]);    
    % Exclude trials with blinks in the saccade-locked time window
    if any(blinkStart>firstSaccade+t_window_saccade(1) & blinkStart<lastSaccade+t_window_saccade(2)) || ...
            any(blinkEnd>firstSaccade+t_window_saccade(1) & blinkEnd<lastSaccade+t_window_saccade(2))
        isBlinky(i) = 1;
    end
    
end

%% Reject eye blink trials
% Remove offending trials from dataset
EEG = pop_rejepoch(EEG,isBlinky,0);
