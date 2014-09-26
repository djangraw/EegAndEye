function EEG = ExcludeTrialsByCutoff(EEG, v_cutoff,t_window_epoch,t_window_saccade,t_discrim,experimentType)

% Takes the current data set in EEGLAB and applies a hard cutoff.  A trial
% on which any electrode exceeds v_cutoff (mV) anywhere within the time 
% window t_window (s) is rejected and removed from the dataset.
%
% EEG = ExcludeTrialsByCutoff(EEG, v_cutoff,t_window_epoch, t_window_saccade,t_discrim,experimentType)
%
% INPUTS:
% -EEG is an eeglab data structure.
% -v_cutoff is the voltage (mV) above which is considered an artifact.
% - t_window_epoch is a 2-element vector specifying times (in ms) of the 
%   epoch in which blinks will affect a trial (and so should be removed). 
%   (default: [-350 0] to avoid baseline)
% - t_window_saccade is a 2-element vector specifying times (in ms) relative
%   to the first and last saccade in an epoch in which blinks will affect a 
%   trial (and so should be removed). (default: [-500 500])
% - t_discrim is a 2-element vector specifying the times (in ms) relative
%   to stimulus onset in which a saccade could be used in discrimination.  
%   Only these saccades should be shielded from high voltages.
% - experimentType is a string indicating the type of experiment
%   ('3DSearch' (default) or 'Squares')
%
% OUTPUT:
% -EEG is the inputted eeglab data structure, with specified trials
% rejected.
%
% Created 8/6/10 by DJ.
% Updated 11/3/10 by DJ - comments
% Updated 11/4/10 by DJ - removed try/catch statement
% Updated 12/10/10 by DJ - switched v_cutoff from 100 to 50
% Updated 2/23/11 by DJ - changed to function. 
% Updated 7/27/11 by DJ - added t_window_saccade input, changed defaults,
% switched from pop_eegthresh to code from RemoveEyeBlinkTrials.
% Updated 8/1/11 by DJ - consider both start and end saccade times, add
%    t_discrim input
% Updated 11/1/11 by DJ - added experimentType input to work with Squares
%    experiments

% Set up
if nargin<2 || isempty(v_cutoff)
    v_cutoff = 50; % max voltage (in mV) allowed in any electrode
end
if nargin<3 
    t_window_epoch = [-200 0]; % time window (in ms) in which to check for this voltage
end
if nargin<4 
    t_window_saccade = [-500 500]; % time window (in ms) in which to check for this voltage
end
if nargin<5
    t_discrim = [0 1000]; % epoch time range (min, max in ms) in which saccades could be used in discrimination.
end     
if nargin<6
    experimentType = '3DSearch';
end
if strcmp(experimentType,'3DSearch')
    GetNumbers;
elseif strcmp(experimentType,'Squares')
    Constants = GetSquaresConstants;
end


%% Apply threshold and reject offending trials
fprintf('---Excluding trials with voltages over %d at t=[%d, %d] or [%d, %d] relative to saccade\n',...
    v_cutoff, t_window_epoch(1),t_window_epoch(2), t_window_saccade(1),t_window_saccade(2))

% sacStartTimes = GetEpochSaccades(EEG,Numbers.SACCADE_START);
% sacEndTimes = GetEpochSaccades(EEG,Numbers.SACCADE_END);

%% Find eye blink trials
isBlinky = zeros(1,numel(EEG.epoch)); % was there a blink in the time window on this trial?
for i=1:numel(EEG.epoch)    
    % 0. Extract necessary info from EEG struct
    % Get events in this epoch
    eventTypes = str2double([EEG.epoch(i).eventtype]);
    eventLatencies = cell2mat([EEG.epoch(i).eventlatency]);
    % Separate out events    
    if strcmp(experimentType,'3DSearch')
        saccadeStart = eventLatencies(eventTypes==Numbers.SACCADE_START);
        saccadeEnd = eventLatencies(eventTypes==Numbers.SACCADE_END);
    elseif strcmp(experimentType,'Squares')
        saccadeStart = eventLatencies(ismember(eventTypes,Constants.SACCADESTART_BASE:Constants.SACCADEEND_BASE));
        saccadeEnd = eventLatencies(ismember(eventTypes,Constants.SACCADEEND_BASE:(2*Constants.SACCADEEND_BASE-Constants.SACCADESTART_BASE)));
    end
    % Get above-threshold times in this epoch
    isAbove = any(EEG.data(:,:,i)>v_cutoff | EEG.data(:,:,i)<-v_cutoff,1);
    timesAbove = EEG.times(isAbove); % times in ms in which voltage of >=1 elec is outside acceptable range
    
    % 1. Exclude trials exceeding threshold in the specified epoch time window
    if  any(timesAbove>t_window_epoch(1) & timesAbove<t_window_epoch(2))
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
    
    % 3. Exclude trials exceeding threshold in the specified saccade-locked time window
    % Get first and last saccades
    firstSaccade = min([okSaccadeStart, okSaccadeEnd]);  
    lastSaccade = max([okSaccadeStart, okSaccadeEnd]);    
    % Exclude trials exceeding threshold in the saccade-locked time window
    if any(timesAbove>firstSaccade+t_window_saccade(1) & timesAbove<lastSaccade+t_window_saccade(2))
        isBlinky(i) = 1;
    end
    
end

%% Reject eye blink trials
% Remove offending trials from dataset
EEG = pop_rejepoch(EEG,isBlinky,0);
EEG = eeg_checkset( EEG );
