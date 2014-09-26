% PlotBionavFigures_03302012.m
%
% This script plots the results from subjects 18, 19, and 20 as calculated
% by RunResultsThroughTag and CalculateImprovement.
% The figures were made for Paul's presentation at the CaN CTA meeting at
% ARL in early April 2012.
%
% Created 3/30/12 by DJ.

% Input results
EEG_pc = [37.5 62.5 25];
TAG_pc = [90 92.1 2.4];
train_pf = [45/600 54/600 63/720]*100;
EEG_pf = [2 3.4 1.2];
TAG_pf = [82.9 87.8 2.4];
pct_dist = [46.8 46.3 46.4];

% Plot precision
figure(1); clf; hold on;
plot([1 1 1;2 2 2],[EEG_pc;TAG_pc],'.-');
plot([1 2], [25 25],'k.--');
% Annotate plot
title('Precision of Predicted Targets');
ylabel('Precision (%)')
axis([0 3 0 100])
set(gca,'xtick',[1 2],'xticklabel',{'EEG','TAG'});
legend('S18','S19','S20','Chance','Location','SouthEast');

% Plot percent found
figure(2); clf; hold on;
plot([1 1 1;2 2 2],[train_pf;TAG_pf],'.-');
plot([1 2], [100/12 25],'k.--');
% Annotate plot
title('Percentage of True Targets Identified');
ylabel('Targets Identified (%)')
axis([0 3 0 100])
set(gca,'xtick',[1 2],'xticklabel',{'Training','TAG'});
legend('S18','S19','S20','Chance','Location','SouthEast');

% Plot search efficiency
figure(3); clf; hold on;
plot([0 0 0; pct_dist],[0 0 0; TAG_pf],'.-');
plot([0 100], [0 100],'k.--');
% Annotate plot
title('Search Efficiency with Full System');
xlabel('% Distance Traveled (normalized)')
ylabel('% Targets Viewed')
axis([0 100 0 100])
legend('S18','S19','S20','View All','Location','SouthEast');