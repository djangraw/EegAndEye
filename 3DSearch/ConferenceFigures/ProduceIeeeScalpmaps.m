% ProduceIeeeScalpmaps
%
% Plots the scalpmaps at given time points for the subjects indicated 
% in the code.  Targets on top row, distractors on the bottom.
% This plot was Figure 3 in DJ's 2011 IEEE/EMBS NE conference submission.
%
% Created 1/14/11 by DJ.
% Updated 1/18/11 by DJ - comments.


% set up
% subjects = [6 2 7]; % all subjects
subjects = [8 10 11 12];
timepoints = [ .150 .450 .750]; % times (s)
color_limits = [-8 8]; % scale of topoplot
figure;
nCols = numel(timepoints);
set(gcf,'Position',[0 200 550 400])

% Load Data
LoadAllEpochs; % alter this program to get stimulus-locked or saccade-locked analysis.

% Calculate ERPs
[erp1 erp2 erpsub time] = pop_comperp( ALLEEG, 1,1:2:(2*numel(subjects)),2:2:(2*numel(subjects)),'addavg','on','subavg','on','diffavg','off','diffstd','off','lowpass',30,'tplotopt',{'ydir',1});
close;


% Plot
for i=1:numel(timepoints)
    % Get time point
    iTime = find(time>=timepoints(i),1);
    % Target plot
    subplot(2,nCols,i);
    topoplot(erp1(:,iTime),ALLEEG(1).chanlocs,'MapLimits',color_limits,'conv','on');
    title(sprintf('Targets, t = %g ms',time(iTime)*1000));  
    % Distractor plot
    subplot(2,nCols,numel(timepoints)+i);
    topoplot(erp2(:,iTime),ALLEEG(1).chanlocs,'MapLimits',color_limits,'conv','on');
    title(sprintf('Distractors, t = %g ms',time(iTime)*1000));  
end
colorbar;
clear timepoints color_limits nRows iTime i