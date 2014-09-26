function PlotLooResults(LOO)

% function PlotLooResults(LOO)
%
% Plots the Az scores from a LOO struct like that created in SaveLooResults
% in a way like that created in MakeLogRegMovie.
%
% Created 1/11/11 by DJ.

%% Plot LOO Results

% fprintf('Std dev = %f\n', std(LOO.Az));
cla; hold on
plot(LOO.time,LOO.Az);
plot(get(gca,'XLim'),[0.5 0.5],'k--');
plot(get(gca,'XLim'),[0.75 0.75],'k:');
ylim([0.3 1]);
plot([0 0],get(gca,'YLim'),'k-');
title(sprintf('%s\nvs. %s',LOO.setname1,LOO.setname2));
xlabel('time (s)');
ylabel('LOO Az');