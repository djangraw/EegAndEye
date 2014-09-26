function GetEyeAverages(ALLEEG,trialType)

% Like AnalyzeEyeEvents, but for multiple subjects.
%
% GetEyeAverages(ALLEEG,trialType)
%
% Plots the times when saccades occurred, on average for all subjects.
% Note that this function should be called after LoadEpochs has put the
% data into EEGLAB.
% This plot was Figure 5 in DJ's 2011 IEEE/EMBS NE conference submission.
%
% INPUTS:
% -ALLEEG is EEGLAB's n-element vector of structs.
% -trialType is an n-element vector in which each dataset is assigned a
%  label (1 for target, 0 for distractor, currently).
%
% Created 1/13/11 by DJ.
% Updated 1/18/11 by DJ - comments.
% Updated 6/2/11 by DJ - switched from SACCADE_TO to SACCADE_END event code

GetNumbers;
tBins = -450:100:950;

% get histograms and plot
for i=1:numel(ALLEEG)
    EEG=ALLEEG(i);
    nTrials(i) = EEG.trials;

    eventTypes = str2double([EEG.epoch(:).eventtype]);
    eventLatencies = cell2mat([EEG.epoch(:).eventlatency]);
    
    % get hist for blinks
%     blinkLatencies{i} = eventLatencies(eventTypes == Numbers.BLINK);
%     yBlink{i} = hist(blinkLatencies{i}/1000,tBins);
%     plot(tBins,yBlink{i}/nTrials(i),colors((i-1)*2+1));
    
    % get hist for saccades
    saccadeLatencies{i} = eventLatencies(eventTypes == Numbers.SACCADE_END);
    ySaccade(i,:) = hist(saccadeLatencies{i},tBins);
%     plot(tBins,ySaccade{i}/nTrials(i),colors((i-1)*2+2));
end


cla; hold on;
colors = 'brgcmy';
trialTypes = unique(trialType);
for i=1:numel(trialTypes)
    isThisType = trialType==trialTypes(i);
    yThisType = sum(ySaccade(isThisType,:),1);
    nThisType = sum(nTrials(isThisType));
    plot(tBins,yThisType/nThisType*100,colors(i));
end

% Annotate plot
plot([0 0],get(gca,'ylim'),'k') % plot zero line
title('Timing of eye events in an epoch')
xlabel('time (ms)')
ylabel('% of trials');
xlim([-200 1000]);
legend('Targets','Distractors')
