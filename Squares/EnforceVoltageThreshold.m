function EEG = EnforceVoltageThreshold(EEG,vThresh,regressor_events,artifact_events,extent_ms,artifact_extent_ms)

% Finds electrodes exceeding a voltage threshold, and interpolates their 
% values in offending epochs.  Marks epochs in which >5 electrodes exceed 
% threshold for rejection.
%
% EEG = EnforceVoltageThreshold(EEG,vThresh,regressor_events,artifact_events,extent_ms,artifact_extent_ms)
%
% INPUTS: 
% -EEG0 is an eeglab data struct with 'FixOn-T' and 'FixOn-D' events
% marking the beginning of each trial and 'BlinkStart' and 'BlinkEnd' 
% events marking the start and end of each blink.
% -tWindow is a 2-element vector marking the start and end time (in s) of
% each epoch relative to the 'TrialStart' event.
% -vThresh is a scalar indicating the threshold (in uV) above which data
% should be considered an artifact.
% -blinkexpand is a 2-element vector indicating the time (in s) outside the 
% blink that should be considered part of the blink (see GetBlinkMask.m).
% For example, when running a GLM in which regressors within 0.5s of a
% blink will be zeroed out, set blinkexpand = [0.5 0.5].
% 
% OUTPUTS: 
% -EEG is an epoched eeglab data struct with data cleaned as described
% above.
% 
% Created 8/8/12 by DJ.
% Updated 3/22/13 by DJ - use asymmetric influence, GetGlmRegressors_v2p0

if nargin<2 || isempty(vThresh)
    vThresh = 75;
end
if nargin<6 || isempty(artifact_extent_ms)
    warning('artifact_extent_ms not specified!');
    artifact_extent_ms = extent_ms;
end


nBadChan_thresh = 5;

dt = 1000/EEG.srate; % in ms
regressor_range = round(extent_ms/dt); % how many samples should each response function extend?
artifact_range = round(artifact_extent_ms/dt); % how many samples should each artifact affect?
t = (1:(EEG.pnts*EEG.trials))*dt; % in ms

event_times = [EEG.event(ismember({EEG.event.type},regressor_events)).latency]*dt;
artifact_times = [EEG.event(ismember({EEG.event.type},artifact_events)).latency]*dt;
% [~,S] = GetGlmRegressors(t,{event_times},artifact_times,Nt); % make sure zeroPadS option is on
[~,S] = GetGlmRegressors_v2p0(t,{event_times},artifact_times,regressor_range,artifact_range); % make sure zeroPadS option is on
isApplicable = any(S,2);
isApplicable = reshape(isApplicable,[EEG.pnts, EEG.trials])';

% Apply voltage Threshold
isBadChan = false(EEG.trials,EEG.nbchan);
for i=1:EEG.trials
    isOver = (abs(EEG.data(:,:,i))>vThresh);
    isOver = isOver & repmat(isApplicable(i,:),EEG.nbchan,1);
    
    % Find offending spots
    for j=1:EEG.nbchan
        if any(isOver(j,:))                        
            isBadChan(i,j) = true;
%             % Find time points where threshold is exceeded
%             ups = find(diff(isOver(j,:))>0);
%             downs = find(diff(isOver(j,:))<0);            
%             if isOver(j,1), ups = [1 ups]; end
%             if isOver(j,end), downs = [downs EEG.pnts]; end
%             % Print results
%             fprintf('Trial %d, Electrode %s: ',i,EEG.chanlocs(j).labels);
%             for k=1:numel(ups)
%                 fprintf('[%.0f %.0f] ',EEG.times(ups(k)), EEG.times(downs(k)));
%             end
%             fprintf('\n');
        end        
    end
    
end


nBadChan = sum(isBadChan,2);
fprintf('%d/%d = %.1f%% trials with >%d bad channels\n',...
    sum(nBadChan>nBadChan_thresh),EEG.trials, ...
    mean(nBadChan>nBadChan_thresh)*100, nBadChan_thresh);
fprintf('%d/%d = %.1f%% trials with 1-%d bad channels\n',...
    sum(nBadChan>=1 & nBadChan<=nBadChan_thresh), EEG.trials, ...
    mean(nBadChan>=1 & nBadChan<=nBadChan_thresh)*100, nBadChan_thresh);

figure;
subplot(2,1,1); cla; hold on;
plot(sum(isBadChan,2),'b.');
plot(get(gca,'xlim'),[5.5 5.5],'k--');
xlabel('trial number')
ylabel('number of bad electrodes')
title(sprintf('Electrodes exceeding %guV in each epoch',vThresh));

subplot(2,1,2); cla;
plot(sum(isBadChan,1),'b.');
% plot(get(gca,'xlim'),[5 5],'k--');
xlabel('electrode number')
ylabel('number of bad trials')
title(sprintf('Epochs exceeding %guV for each electrode',vThresh));


% interpolate electrodes on offending trials
fprintf('Interpolating channels in epochs where they exceed %guV...\n',vThresh)
for i=1:EEG.trials    
    EEG = InterpChan(EEG,EEG.chanlocs(isBadChan(i,:)),i); 
end

fprintf('Removing epochs where >%d electrodes exceed threshold...\n',nBadChan_thresh)
% Mark bad trials for rejection
isBadTrial = (nBadChan>nBadChan_thresh);
if ~isfield(EEG.etc,'rejectepoch')
    EEG.etc.rejectepoch = isBadTrial;
else
    EEG.etc.rejectepoch = EEG.etc.rejectepoch | isBadTrial;
end
% Reject bad trials
% EEG = pop_rejepoch(EEG,iBadTrials,0);