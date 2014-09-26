function [n,xhist] = FindTimeHistOfEvents(R,event_types,xhist);

% [n,xhist] = FindTimeHistOfEvents(R,event_types,xhist);
%
% Multi-subject saccade time histogram.
%
% Created 1/28/13 by DJ for one-time use.

if nargin<3 || isempty(xhist)
    xhist = 25:50:4000;
end
nEvents = numel(event_types);

% clear n
for iSubj = 1:numel(R) % subject       
    for i=1:nEvents % square
        % get
        tb = FindTimeBetweenEvents(R(iSubj).EEG,'TrialStart-D',event_types{i});
        tb2 = FindTimeBetweenEvents(R(iSubj).EEG,'TrialStart-T',event_types{i});
        % combine
        tb(~isnan(tb2)) = tb2(~isnan(tb2));
        % histogram
        n(i,:,iSubj)=hist(tb,xhist)/sum(~isnan(tb));
        % store mean
        meantime(i,iSubj) = nanmean(tb);
    end
end
% average across subjects
n_avg = mean(n,3);

%% Plot
plot(xhist,n_avg'*100)
xlabel('time from Trial Start (ms)');
ylabel('% of trials');
title(sprintf('Histogram of event times in trial\n(mean across %d subjects)',numel(R)))
legend(event_types);
hold on
PlotVerticalLines(mean(meantime,2),'k--');
    


% %% LSaccade vs. RSaccade
% xhist = 25:50:4000;
% clear n
% for iSubj = 1:9 % subject       
%     % get LSaccade
%     tb = FindTimeBetweenEvents(R(iSubj).EEG,'TrialStart-D','LSaccade');
%     tb2 = FindTimeBetweenEvents(R(iSubj).EEG,'TrialStart-T','LSaccade');
%     % combine
%     tb(~isnan(tb2)) = tb2(~isnan(tb2));
%     % histogram
%     n(1,:,iSubj)=hist(tb,xhist)/sum(~isnan(tb));
%     
%     % get RSaccade
%     tb = FindTimeBetweenEvents(R(iSubj).EEG,'TrialStart-D','RSaccade');
%     tb2 = FindTimeBetweenEvents(R(iSubj).EEG,'TrialStart-T','RSaccade');
%     % combine
%     tb(~isnan(tb2)) = tb2(~isnan(tb2));
%     % histogram
%     n(2,:,iSubj)=hist(tb,xhist)/sum(~isnan(tb));
% end
% % average across subjects
% n_avg = mean(n,3);
% 
% %% Plot
% plot(xhist,n_avg'*100)
% xlabel('time from Trial Start (ms)');
% ylabel('% of trials');
% title(sprintf('Histogram of saccade times in trial\n(mean across 9 subjects)'))
% legend('LSaccade','RSaccade');

