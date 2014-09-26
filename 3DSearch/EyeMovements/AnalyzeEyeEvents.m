% AnalyzeEyeEvents
%
% Takes two of the datasets output by EpochEeglabFile and outputs info 
% about the blinks and saccades associated with those epochs.  
% - First plots a histogram of times at which blinks and saccades occur
% relative to the anchoring events in each epoch.
% - Then displays a few statistics that speak to whether or not the blinks
% and saccades are significantly correlated with one epoch or the other in
% a way that might affect results in other analyses.
%
% Created 11/4/10 by DJ.
% Updated 11/12/10 by DJ - comments
% Updated 12/16/10 by DJ - divide by nTrials, not nSaccades
% Updated 6/2/11 by DJ - switched from SACCADE_TO to SACCADE_END event code


%% Set up
whichDatasets = [4 5]; % [2 3] for target appear and distractor appear (see EpochEeglabFile)
colors = 'brcmgy';

%% Epoch the data
filename = sprintf('3DS-%d-all-filtered.set',subject);
GetNumbers;
epoch_times = [-500 1000]; % time, in ms, that each epoch should extend around the anchoring event
baseline_times = [-200 0]; % time, in ms, relative to the anchoring event, that should be considered baseline.
event_numbers = [Numbers.ENTERS+Numbers.TARGET, Numbers.ENTERS+Numbers.DISTRACTOR,...
    Numbers.SACCADE_TO+Numbers.TARGET, Numbers.SACCADE_TO+Numbers.DISTRACTOR];
event_names = {'targapp','distapp','targsac','distsac'}; % used to name new datasets
refChans = 'None';

EpochEeglabFile(filename,epoch_times,baseline_times,event_numbers,event_names,refChans);

%% Make histograms
% Set up histogram
GetNumbers;
figure; hold on;

% get histograms and plot
for i=1:numel(whichDatasets)
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',whichDatasets(i),'study',0); 
    nTrials(i) = EEG.trials;
    tBins = linspace(EEG.xmin,EEG.xmax,20);
    eventTypes = str2double([EEG.epoch(:).eventtype]);
    eventLatencies = cell2mat([EEG.epoch(:).eventlatency]);
    
    % get hist for blinks
    blinkLatencies{i} = eventLatencies(eventTypes == Numbers.BLINK);
    yBlink{i} = hist(blinkLatencies{i}/1000,tBins);
    plot(tBins,yBlink{i}/nTrials(i),colors((i-1)*2+1));
    
    % get hist for saccades
    saccadeLatencies{i} = eventLatencies(eventTypes == Numbers.SACCADE_END);
    ySaccade{i} = hist(saccadeLatencies{i}/1000,tBins);
    plot(tBins,ySaccade{i}/nTrials(i),colors((i-1)*2+2));
end

% Annotate plot
title('Timing of eye events in an epoch')
xlabel('time from locked event (s)')
ylabel('percentage of trials');
% legend(sprintf('%s - %d blinks',ALLEEG(whichDatasets(1)).setname, numel(blinkLatencies{1})), ...
%     sprintf('%s - %d blinks',ALLEEG(whichDatasets(2)).setname, numel(blinkLatencies{2})))
legend(sprintf('%s - %d blinks',ALLEEG(whichDatasets(1)).setname, numel(blinkLatencies{1})), ...
    sprintf('%s - %d saccades',ALLEEG(whichDatasets(1)).setname, numel(saccadeLatencies{1})), ...
    sprintf('%s - %d blinks',ALLEEG(whichDatasets(2)).setname, numel(blinkLatencies{2})), ...
    sprintf('%s - %d saccades',ALLEEG(whichDatasets(2)).setname, numel(saccadeLatencies{2})));

%% Run statistical tests
disp('...Running statistics...');
% Are there more blinks associated with targets than distractors?
expectedRatio = numel(ALLEEG(whichDatasets(1)).epoch)/( numel(ALLEEG(whichDatasets(1)).epoch) ...
    + numel(ALLEEG(whichDatasets(2)).epoch) );
blinkRatio = numel(blinkLatencies{1})/( numel(blinkLatencies{1}) + numel(blinkLatencies{2}) );
saccadeRatio = numel(saccadeLatencies{1})/( numel(saccadeLatencies{1}) + numel(saccadeLatencies{2}));
fprintf('TR = nNearTarget / (nNearTarget + nNearDistractor) should be %0.3f:\n', expectedRatio);
fprintf('  blinks: TR = %0.3f\n',blinkRatio);
fprintf('  saccades: TR = %0.3f\n',saccadeRatio);

% Is the distribution of each different from a uniform distribution?
expectedMedian = mean(tBins);
pTargBlink = signrank(blinkLatencies{1},expectedMedian);
pDistBlink = signrank(blinkLatencies{2},expectedMedian);
pTargSaccade = signrank(saccadeLatencies{1},expectedMedian);
pDistSaccade = signrank(saccadeLatencies{2},expectedMedian);
fprintf('Median latency is %0.1f s:\n', expectedMedian);
fprintf('  target blinks: p = %0.3f\n',pTargBlink);
fprintf('  distractor blinks: p = %0.3f\n',pDistBlink);
fprintf('  target saccades: p = %0.3f\n',pTargSaccade);
fprintf('  distractor saccades: p = %0.3f\n',pDistSaccade);

% Are the distributions different from each other?
pBlinkDiff = ranksum(blinkLatencies{1},blinkLatencies{2});
pSaccadeDiff = ranksum(saccadeLatencies{1},saccadeLatencies{2});
fprintf('distributions have same medians:\n')
fprintf('   blinks: p = %0.3f\n',pBlinkDiff);
fprintf('   saccades: p = %0.3f\n',pSaccadeDiff);

%% Clean up
clear whichDatasets colors Numbers tBins eventTypes eventLatencies blinkLatencies saccadeLatencies *Ratio p*Blink* p*Saccade* %yBlink ySaccade nTrials