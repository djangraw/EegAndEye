% GetTimingOffsets_script.m
%
% Created 7/21/14 by DJ.

experiment = 'sq';
subject = 17;
suffix = '-all-40Hz-fs100';
EEG = pop_loadset(sprintf('%s-%d%s.set',experiment,subject,suffix));
tWin = [-.1 .1];
EEG = pop_epoch(EEG,{'BS'},tWin);
template = mean(EEG.data,3);
MakeTopoMovie(double(template),EEG.times,EEG.chanlocs);

%%
experiment = 'sf3';
subjects = 1:12;

% experiment = 'sf';
% subjects = [1:10 12:13];

% experiment = 'sq';
% subjects = [9:11, 13:15, 17:19, 20:27];

clear erp;
for i=1:numel(subjects)
    subject = subjects(i);
    
    EEG = pop_loadset(sprintf('%s-%d%s.set',experiment,subject,suffix));
    EEG = pop_epoch(EEG,{'BS'},[-.5 .5]);
    erp(:,:,i) = mean(EEG.data,3);
    
end

%% Use Woody algorithm to line up trials 
times = EEG.times;
tWin = [-.1 .1];
windowSize = size(template,2);
t = times((1:end-windowSize)+round(windowSize/2))-mean(tWin)*1000;
matchStrength = UpdateTemplateMatchStrength(erp,template);

%% Plot results

clf; hold on; 
imagesc(t,1:numel(subjects),matchStrength);

[~,iMax] = max(matchStrength,[],2);

for i=1:numel(subjects)
    plot(t(iMax(i)),i,'k.');
end
axis([t(1) t(end) .5 numel(subjects)+.5])
    
