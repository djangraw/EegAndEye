function EEG = RemoveIncorrectTrials(EEG)

% Reject target epochs in which a button was not pressed and distractor
% epochs in which a button was pressed.
%
% EEG = RemoveIncorrectTrials(EEG)
%
% INPUTS:
% -EEG is a standard EEGLAB data struct from a 3DSearch experiment.
%
% OUTPUTS:
% -EEG is the same as the input but with offending trials removed.
% 
% Created 8/4/11 by DJ.

% Display header
disp('---Removing Incorrect Trials...');

% extract event info
eventTypes = [str2double({EEG.event(:).type})];
eventLatencies = [EEG.event(:).latency] * 1000/EEG.srate;
eventEpochs = [EEG.event(:).epoch];

% extract epoch info
GetNumbers;
isTargetEpoch = zeros(1,EEG.trials);
isTargetEpoch(eventEpochs(eventTypes==Numbers.TARGET+Numbers.ENTERS)) = 1;
isButtonEpoch = zeros(1,EEG.trials);
isButtonEpoch(eventEpochs(eventTypes==Numbers.BUTTON)) = 1;

% find bad trials
badTargetTrials = find(isTargetEpoch & ~isButtonEpoch);
badDistractorTrials = find(~isTargetEpoch & isButtonEpoch);

% inform user
if sum(isTargetEpoch)>0
    if numel(badTargetTrials)==sum(isTargetEpoch) 
        error('All target trials marked as incorrect!  Make sure button events were logged correctly.');
    end
    fprintf('Excluding %d/%d target trials with no button press...\n',numel(badTargetTrials),sum(isTargetEpoch));       
end
if sum(~isTargetEpoch)>0
    if numel(badDistractorTrials)==sum(~isTargetEpoch) 
        error('All distractor trials marked as incorrect!  Make sure button events were logged correctly.');
    end
    fprintf('Excluding %d/%d distractor trials with button presses...\n',numel(badDistractorTrials),sum(~isTargetEpoch));
end

% remove trials
EEG = pop_rejepoch(EEG,[badTargetTrials badDistractorTrials],0);