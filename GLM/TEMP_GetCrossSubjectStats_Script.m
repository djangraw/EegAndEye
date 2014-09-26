% Created 7/13/12 by DJ.
% Updated 7/18/12 by DJ - all one figure.
% Updated 3/19/14 by DJ - tResponse in cells.

%% CALCULATE
nSubjects = size(timecourse{1},1);
nT = size(timecourse{1},2);
difference = nan(nComponents,nT,nTrialTypes,nTrialTypes); % sign of significant differences between classes
disp('Getting p values...');
hwait = waitbar(0,'Getting p values...');
for j=1:nT %nTimepoints
    waitbar(j/nT,hwait);
    for i=1:nComponents   
        X = nan(nSubjects,nTrialTypes);
        for k=1:nTrialTypes
            X(:,k) = timecourse{i,k}(:,j);
        end
        [~,~,stats] = anova1(X,S(1).regressor_events{S(1).iLevel},'off');
        comparison = multcompare(stats,'alpha',0.05,'display','off','ctype','tukey-kramer');
        sig_pos = find(comparison(:,3)>0);
        sig_neg = find(comparison(:,5)<0);
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
close(hwait);
disp('Done!')

%% PLOT

colors = {'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r--' 'g--' 'b--'};
if iscell(S(1).tResponse)
    tResponse = S(1).tResponse{S(1).iLevel};
else
    tResponse = S(1).tResponse;
end
figure(55); clf;
for j=1:nTrialTypes
%     figure(55+j); clf;
    for i = 1:nComponents
        subplot(nComponents,nTrialTypes,(i-1)*nTrialTypes+j)
        cla; hold on;
        for k=1:nTrialTypes
            plot(tResponse,difference(i,:,j,k)*k,colors{k});
        end
        axis([tResponse(1), tResponse(end), -nTrialTypes-1, nTrialTypes+1]);
        set(gca,'xgrid','on')
        plot(get(gca,'xlim'),[0 0],'k--');
        ylabel(sprintf('Component %d',i));
    end
    xlabel('time (ms)');
    subplot(nComponents,nTrialTypes,j);
    title(sprintf('%s vs. all others (One-Way Anova with multcompare, p=0.05)',S(1).regressor_events{S(1).iLevel}{j}));
    subplot(nComponents,nTrialTypes,1);
    legend(S(1).regressor_events{S(1).iLevel},'Location','NorthWest');
end
