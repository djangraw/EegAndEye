% PlotEyeAndEegForSquaresTrial
% Created 2/16/12 for one-time use by DJ.

trialNum = 6;
elRange = [y(1).trial.fix_time(trialNum), y(1).trial.circle_time(trialNum)];
PlotEyeSamplesAndEvents(eyepos,29478,y(1).saccade.start_time,'saccadeStart',...
    elRange);

% hEOG approximated by AF7 (chan 5) minus AF8 (chan 13)
eegRange = round(EyelinkToEegTimes(elRange,y(1)));
plot((0:eegRange(2)-eegRange(1))/EEG.srate*1000, ...
    EEG.data(5,eegRange(1):eegRange(2))-EEG.data(13,eegRange(1):eegRange(2)),'r')