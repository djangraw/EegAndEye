function [tPup, dPup] = PlotPupilErps(data, baselineWin)

% [tPup, dPup] = PlotPupilErps(data);
%
% INPUTS:
% -data is the data field of the output struct of GetEyeErps.
% -baselineWin is a 2-element vector indicating the start and end times of
% the baseline time (typically just before an object appears).
%
% OUTPUTS:
% -tPup is a matrix of the pupil diameter over time for target trials.
% -dPup is a matrix of the pupil diameter over time for distractor trials.
%
% Created 6/20/11 by DJ.


% Handle defaults
if nargin<2
    baselineWin = [-200 0];
end

% Extract info
tPup = data.targetPupEpochs;
dPup = data.distractorPupEpochs;
times = data.epochTimes;
% Find baseline
isBaseline = times>=baselineWin(1) & times<=baselineWin(2);

% Subtract baseline
for i=1:size(tPup,1)
    % subtract baseline    
    trialBaseline = mean(tPup(i,isBaseline));
    tPup(i,:) = tPup(i,:) - trialBaseline;
end

for i=1:size(dPup,1)
    % subtract baseline    
    trialBaseline = mean(dPup(i,isBaseline));
    dPup(i,:) = dPup(i,:) - trialBaseline;
end

% Get means and stddevs
tMean = nanmean(tPup);
dMean = nanmean(dPup);
tStd = nanstd(tPup)/sqrt(size(tPup,1));
dStd = nanstd(dPup)/sqrt(size(dPup,1));

% Plot
cla; hold on;
plot(times,tMean,'b','linewidth',2)
plot(times,dMean,'r','linewidth',2)
JackKnife(times,tMean,tStd,'b','b');
JackKnife(times,dMean,dStd,'r','r');

% Annotate
xlabel('time from stimulus onset (ms)')
ylabel('Pupil diameter (a.u.)')
title('Targets vs. Distractors')
legend('targets +/- stderr','distractors +/- stderr','Location','SouthWest')
