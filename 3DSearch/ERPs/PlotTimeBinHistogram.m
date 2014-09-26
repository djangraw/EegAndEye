function PlotTimeBinHistogram(allSigTimes, tBins, divideBy)

% PlotTimeBinHistogram(allSigTimes,tBins, divideBy)
%
% First use LoadEpochs to get results for all subjects.
% Then use PlotGroupErps to get the times when target and distractor trials
% are significantly separated (for each electrode).
% Then run this program.
%
% Created 1/11/11 by DJ.
% Updated 1/18/11 by DJ - comments

% Set up
if nargin<2 || isempty(tBins)
    tBins = -450:100:950;
end
if nargin<3 || isempty(divideBy)
    divideBy = 1;
end

% Get histogram
n = hist(allSigTimes,tBins); % number of electrodes significant in this time bin
pctSig = n/divideBy*100; % avg percent of electrodes significant in this time bin

% Plot
bar(tBins,pctSig)
xlim([-200 1000])
set(gca,'xtick',-200:200:1000)
xlabel('time (ms)');
ylabel('% of electrodes')