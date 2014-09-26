function difference = GetCrossSubjectStats(timecourse,tResponse,regressor_events,Cmap)

% Find and plot the time points at which any two experimental conditions
% evoked a significantly different timecourse of activity.
%
% difference = GetCrossSubjectStats(timecourse,tResponse,regressor_events)
%
% Uses anova1 and multcompare, alpha = 0.05 and tukey-kramer
% multiple-comparisons correction, for statistical tests.
%
% INPUTS:
% - timecourse is an nComponents x nTrialTypes cell matrix, each of which
% contains an nSubjects x nTimePoints matrix of doubles.
% timecourse{i,j}(k,l) is the mean activity in component i evoked by event
% type j for subject k at time point l.  It is an output of
% GetGroupSvdResults.m.
% - tResponse is a nTimePoints-element array of latencies (in ms) relative
% to the locking event.
% - regressor_events is a nTrialTypes cell matrix of strings indicating the
% event types you are comparing.
%
% OUTPUTS:
% - difference is a matrix of size (nComponents,nTimePoints,nTrialTypes,
% nTrialTypes).  difference(i,j,k,l)=1 if, in component i at time j, event
% l evoked significantly more positive activity than event k.  =-1 if more
% negative, and =NaN otherwise.
%
% Created 7/13/12 by DJ.
% Updated 7/18/12 by DJ - Made into function, plots all on one figure.
% Updated 10/2/13 by DJ - Added Cmap input, floating legend

%% CALCULATE
% Get constants
nComponents = size(timecourse,1);
nTrialTypes = size(timecourse,2);
nSubjects = size(timecourse{1},1);
nT = size(timecourse{1},2);

if nargin<4 || isempty(Cmap)
    Cmap = distinguishable_colors(nTrialTypes,{'w','k'}); % don't allow black or white
end


% Get p values using anova1 and multcompare
difference = nan(nComponents,nT,nTrialTypes,nTrialTypes); % sign of significant differences between classes
disp('Getting p values...');
hwait = waitbar(0,'Getting p values...');
for j=1:nT %nTimepoints
    waitbar(j/nT,hwait);
    for i=1:nComponents 
        % Perform statistical test
        X = nan(nSubjects,nTrialTypes); % reformat for anova1
        for k=1:nTrialTypes
            X(:,k) = timecourse{i,k}(:,j);
        end
        [~,~,stats] = anova1(X,regressor_events,'off');
        comparison = multcompare(stats,'alpha',0.05,'display','off','ctype','tukey-kramer'); % correct for multiple comparisons (across conditions, not time points!)
        
        % Find points with significant differences (in positive and negative direction separately)
        sig_pos = find(comparison(:,3)>0);
        sig_neg = find(comparison(:,5)<0);
        % Register these points in the difference matrix
        for k=1:numel(sig_pos)
            difference(i,j,comparison(sig_pos(k),1),comparison(sig_pos(k),2)) = -1;
            difference(i,j,comparison(sig_pos(k),2),comparison(sig_pos(k),1)) = 1;
        end
        for k=1:numel(sig_neg)
            difference(i,j,comparison(sig_neg(k),1),comparison(sig_neg(k),2)) = 1;
            difference(i,j,comparison(sig_neg(k),2),comparison(sig_neg(k),1)) = -1;
        end    
    end
end
% Finish up
close(hwait);
disp('Done!')

%% PLOT

% colors = {'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r--' 'g--' 'b--'};

% figure(55); 
clf;
for j=1:nTrialTypes    
    for i=1:nComponents
        % Plot lines where differences are significant
        subplot(nComponents,nTrialTypes,(i-1)*nTrialTypes+j)
        cla; hold on;
        for k=1:nTrialTypes
            plot(tResponse,difference(i,:,j,k)*k,'color',Cmap(k,:),'linewidth',2);
        end
        % Annotate plot
        plot([tResponse(1) tResponse(end)],[0 0],'color',Cmap(j,:),'linewidth',2);
        axis([tResponse(1), tResponse(end), -nTrialTypes-1, nTrialTypes+1]);
        set(gca,'xgrid','on')
        plot(get(gca,'xlim'),[0 0],'k:','linewidth',2);
        ylabel(sprintf('Component %d',i));
    end
    % Annotate plots futher    
    xlabel('time (ms)'); % only on bottom row
    subplot(nComponents,nTrialTypes,j);
    title(sprintf('%s vs. all others\n(One-Way Anova with multcompare, p=0.05)',regressor_events{j})); % only on top row
    subplot(nComponents,nTrialTypes,1);
%     legend(regressor_events,'Location','NorthWest'); % only on top-left plot
end
MakeLegend(Cmap,regressor_events,2);
