% SaveEpochs.m
%
% Variable 'subject' must be defined.
%
% Created 12/21/10 by DJ.
% Updated 3/15/11 by DJ - EpochEeglabFile and RemoveEyeBlinkTrials are now functions.
% Updated 7/27/11 by DJ - work with new version of RemoveEyeBlinkTrials and
% ExcludeTrialsByCutoff.

%% EPOCH
% subject = 6;
% EpochEeglabFile('3DS-2-all-filtered');
GetNumbers;
[ALLEEG, EEG, CURRENTSET] = EpochEeglabFile(sprintf('3DS-%d-all-filtered-noduds.set',subject),[-500 1000],[-200 0],Numbers.ENTERS+[Numbers.TARGET, Numbers.DISTRACTOR],{'targapp','distapp'},'None');

%% CLEAN
for j=2:3 % the sets that will be saved later
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',j,'study',0); 
    EEG = RemoveEyeBlinkTrials(EEG,Numbers.BLINK,[-350 1000],[]);    
%         PlotExclusionHistogram(EEG.data,EEG.times,[0 1000]);
%         EEG = ExcludeTrialsByCutoff(EEG,100,[-350 1000],[]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end

%% SAVE
% Save results
% pop_saveset(ALLEEG(2), 'filename',sprintf('%d-TargetSac.set',subject));
% pop_saveset(ALLEEG(3), 'filename',sprintf('%d-DistractorSac.set',subject));
pop_saveset(ALLEEG(2), 'filename',sprintf('%d-TargetApp.set',subject));
pop_saveset(ALLEEG(3), 'filename',sprintf('%d-DistractorApp.set',subject));
eeglab redraw;
