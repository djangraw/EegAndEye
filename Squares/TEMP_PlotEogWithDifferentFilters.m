% TEMP_plotEogWithDifferentFilters
%
% Created 6/28/12 by DJ for one-time use.


% EEG13 is the unfiltered data
% EEG13_noeog is the data with EOG regressed out
% EEG13_filtered is the data with EOG regressed out, filtered from 1-50Hz
% EEG0p1 is the data filered from 0.1-100Hz
% EEG is the data filtered from 1-50Hz

%% Load event data
fixEvents = [EEG.event(ismember({EEG.event.type},{'FixOn-D','FixOn-T'})).latency]/EEG.srate;
t = (1:size(samples.eyepos,1))/1000 + (EyelinkToEegTimes(4757158,x)/x.eeg.eventsamplerate);

%% Plot
iElec = find(strcmp('F8',{EEG13(1).chanlocs.labels}));

clf
subplot(3,1,1)
cla; box on; hold on;
title(sprintf('Electrode %s',EEG.chanlocs(iElec).labels))
plot((1:EEG13(1).pnts)/EEG13(1).srate, EEG13(1).data(iElec,:))
plot(t,samples.eyepos(:,1)/10-200,'g')
ylim([-200 200])
PlotVerticalLines(fixEvents,'r')
plot(get(gca,'xlim'),[0 0],'k')
ylabel('Unfiltered')
legend('EEG','EyePos','TrialStart','Location','NorthWest')

subplot(3,1,2)
cla; box on; hold on;
% plot((1:EEG0p1.pnts)/EEG0p1.srate, EEG0p1.data(23,:))
plot((1:EEG13_noeog.pnts)/EEG13_noeog.srate, EEG13_noeog.data(iElec,:))
plot(t,samples.eyepos(:,1)/10-200,'g')
ylim([-200 200])
PlotVerticalLines(fixEvents,'r')
plot(get(gca,'xlim'),[0 0],'k')
% ylabel('0.1-100Hz')
ylabel('EogRemoved')

subplot(3,1,3)
cla; box on; hold on;
% plot((1:EEG.pnts)/EEG.srate, EEG.data(23,:))
plot((1:EEG13_filtered.pnts)/EEG13_filtered.srate, EEG13_filtered.data(iElec,:))
% plot((1:EEG13_noblinks.pnts)/EEG13_noblinks.srate, EEG13_noblinks.data(iElec,:))
plot(t,samples.eyepos(:,1)/10-200,'g')
ylim([-200 200])
PlotVerticalLines(fixEvents,'r')
plot(get(gca,'xlim'),[0 0],'k')
ylabel('1-50Hz')
xlabel('time (s)')