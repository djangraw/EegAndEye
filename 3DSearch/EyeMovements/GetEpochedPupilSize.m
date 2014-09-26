function [ps_epoch, t_epoch, isTargetEpoch] = GetEpochedPupilSize(y,ps,EEG,epochRange)

% Gets matrix of epoched pupil size values.
%
% [ps_epoch, t_epoch, isTargetEpoch] = GetEpochedPupilSize(y,ps,EEG,epochRange)
% [ps_epoch, t_epoch, isTargetEpoch] = GetEpochedPupilSize(y,ps,EEG,offset)
%
% INPUTS:
% 
% OUTPUTS:
%
% Created 7/8/13 by DJ.

% Handle inputs
if nargin<4 || isempty(epochRange)
    epochRange = [EEG.xmin*1000, EEG.xmax*1000];
elseif numel(epochRange)==1
    offset = epochRange;
    epochRange = [EEG.xmin*1000, EEG.xmax*1000] - offset;
end
t_epoch = epochRange(1):epochRange(2); % assume 1000Hz sampling rate

Numbers = GetNumbers;
target_event_str = num2str(Numbers.SACCADE_TO + Numbers.TARGET);

%% Get initial time of each epoch
zeroTime = nan(1,EEG.trials);
zeroSession = nan(1,EEG.trials);
isTargetEpoch = nan(1,EEG.trials);
urevent_isboundary = strcmp('boundary',{EEG.urevent.type});
urevent_session = 1+cumsum(urevent_isboundary);
for i=1:EEG.trials
    iZero = find([EEG.epoch(i).eventlatency{:}]==0,1);
    zeroTime(i) = EEG.epoch(i).eventinit_time{iZero}; % in seconds
    isTargetEpoch(i) = any(strcmp(target_event_str,EEG.epoch(i).eventtype));
    zeroSession(i) = urevent_session(EEG.epoch(i).eventurevent{iZero});
end
    
%% Get pupil epochs
nSessions = numel(y);
ps_epoch = [];
for i=1:nSessions
    x = y(i);
    pupilsize = ps{i}';    
    tPupil = (1:length(pupilsize)) + x.eyelink.record_time - 1;
    
    % get eyelink time zero
    zeroTime_pupil = round(EyelinkToEegTimes(zeroTime(zeroSession==i),x.eeg.sync_events(:,1)/x.eeg.eventsamplerate,x.eyelink.sync_events(:,1)));
    for j=1:numel(zeroTime_pupil)
        iEpochStart = find(tPupil==zeroTime_pupil(j)+epochRange(1));
        iEpochEnd = find(tPupil==zeroTime_pupil(j)+epochRange(2));
        % Pad with NaN's if necessary
        if isempty(iEpochStart) && isempty(iEpochEnd)
            fprintf('No Pupil info! Session, trial: %d, %d\n',i,j);
            ps_epoch = [ps_epoch; nan(1,length(t_epoch))];
        elseif isempty(iEpochStart)
            fprintf('Padding start of Session, trial: %d, %d\n',i,j);            
            known_data = pupilsize(1:iEpochEnd);
            padded_epoch = nan(1,length(t_epoch));            
            padded_epoch(end-length(known_data)+1:end) = known_data;
            ps_epoch = [ps_epoch; padded_epoch];
        elseif isempty(iEpochEnd)
            fprintf('Padding end of Session, trial: %d, %d\n',i,j);
            known_data = pupilsize(iEpochStart:end);
            padded_epoch = nan(1,length(t_epoch));            
            padded_epoch(1:length(known_data)) = known_data;
            ps_epoch = [ps_epoch; padded_epoch];
        else
            ps_epoch = [ps_epoch; pupilsize(iEpochStart:iEpochEnd)];
        end
    end
        
end



