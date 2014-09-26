function [avgPlus, avgMinus, separatedTimes] = PlotChannelErp(dataPlus, dataMinus, time, channel)

% Plots the ERP for a single channel, including the difference and standard error bars.
% Results will be filtered with a 30-Hz LPF (Butterworth 5th order) before plotting.
%
% [avgPlus, avgMinus, separatedTimes] = PlotChannelErp(dataPlus, dataMinus, time)
%
% INPUTS:
% dataPlus is an mxn matrix of data (in uV), where each row is a trial.
% dataMinus is an mxn matrix of data (in uV), where each row is a trial.
% time is an n-element vector of time values for the data samples (in ms).
% channel is a string or number indicating the channel being plotted.
%
% OUTPUTS:
% -avgPlus and avgMinus are the means of the given data (with filters
% applied).
% -separatedTimes is a vector of time points at which the two datasets were
% significantly separated (green crosses on plot).
%
% Created 12/20/10 by DJ.
% Updated 1/7/11 by DJ - added outputs
% Updated 1/14/11 by DJ - no longer makes a new figure automatically

%% SET UP
% Handle inputs
if nargin<4
    channel = '';
end
nPlus = size(dataPlus,1);
nMinus = size(dataMinus,1);

% Convert to doubles
dataPlus = double(dataPlus);
dataMinus = double(dataMinus);
% dataDiff = dataPlus-dataMinus;

% Calculate averages
avgPlus = mean(dataPlus,1);
avgMinus = mean(dataMinus,1);
% avgDiff = mean(dataDiff,1);
avgDiff = avgPlus - avgMinus;

% Calculate standard error bars
stdPlus = std(dataPlus,1)/sqrt(nPlus);
stdMinus = std(dataMinus,1)/sqrt(nMinus);
% stdDiff = std(dataDiff,1)/sqrt(nRows);

%% FILTER
filter_cutoff = 30;
filter_order = 5;
fs = 1000/(time(2)-time(1));
[b, a] = butter(filter_order,filter_cutoff/fs*2);
% Filter averages
avgPlus = filtfilt(b,a,avgPlus);
avgMinus = filtfilt(b,a,avgMinus);
avgDiff = filtfilt(b,a,avgDiff);

% % Filter standard error bars
% stdPlus = filtfilt(a,b,stdPlus);
% stdMinus = filtfilt(a,b,stdMinus);
% stdDiff = filtfilt(a,b,stdDiff);

%% STATISTICS
alpha = 0.05;
for i=1:size(dataPlus,2)
    [p(i),h(i)] = ranksum(dataPlus(:,i),dataMinus(:,i),'alpha',alpha); % test hypothesis that median is zero
end
separatedTimes = time(h);
fprintf('%d significant epochs at alpha=%g\n',sum(h),alpha);

%% PLOT
% Format axes
cla;
xlim([-200 1000]);
ylim([-10 20]);
box on; hold on;

% Plot Averages and Standard Error Bars
JackKnife(time,avgPlus,stdPlus,'b','b');
JackKnife(time,avgMinus,stdMinus,'r','r');
plot(time,avgDiff,'k','linewidth',2); % without error bars
% JackKnife(time,avgDiff,stdDiff,'k','k'); % with error bars

% Plot significant times
ylimits = get(gca,'ylim');
plot(separatedTimes,repmat(ylimits(2)-1,1,numel(separatedTimes)),'g+');

% Plot axes
plot(get(gca,'xlim'),[0 0],'k');
plot([0 0],get(gca,'ylim'),'k');

% Annotate Plot
xlabel('time (ms)');
ylabel('voltage (uV)');
title(['Channel ' num2str(channel) ' ERP']);
MakeLegend({'b-' 'r-'},{sprintf('Targets (n=%d)',nPlus) sprintf('Distractors (n=%d)',nMinus)},[2 2]);



