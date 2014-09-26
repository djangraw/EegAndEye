% GetMeanSpectra.m
%
% Created for one-time use on 12/16/11

% i = 4; % subj number
chanName = 'FZ';
% freqs  = [0 4 8 12 24 48];
% 
% 
% [S0{i},F0{i},T0{i},P0{i}] = PlotSpectrogram(ALLEEG(2), chanName, freqs);
% [S1{i},F1{i},T1{i},P1{i}] = PlotSpectrogram(ALLEEG(3), chanName, freqs);
% [SI{i},FI{i},TI{i},PI{i}] = PlotSpectrogram(ALLEEG(4), chanName, freqs);
% [SC{i},FC{i},TC{i},PC{i}] = PlotSpectrogram(ALLEEG(5), chanName, freqs);
% 

chanLabels = {EEG.chanlocs(:).labels}; % the names of all the channels
chan = find(strcmp(chanName,chanLabels));

%%

T = T0{1};
F = F0{1};
P = mean(cat(3,PC{:}),3)./mean(cat(3,PI{:}),3);

timeStart = EEG.times(1)/1000;
figure; surf(T + timeStart, F, 10*log10(P),'edgecolor','none'); axis tight; 
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');
title({sprintf('Power (in dB) for channel %d (%s)', chan, EEG.chanlocs(chan).labels); EEG.setname});
colorbar;
ylim([min(F) max(F)]);
caxis([-20 20]);