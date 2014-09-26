function PlotResultsAsErps(group_RF,group_Z,chanlocs,tResponse,legendstr,Cmap,chansToPlot,titlestr)

% PlotResultsAsErps(group_RF,group_Z,tResponse,chanlocs,legendstr,Cmap,chansToPlot,titlestr)
%
% Created 9/16/14 by DJ based on GetContrastResults_fast.
% Updated 9/18/14 by DJ - allow any size chansToPlot

%% Set up
if ~exist('Cmap','var') || isempty(Cmap)
    Cmap = distinguishable_colors(size(group_Z,3),{'w','k'});
end
if iscell(Cmap)
    Cmap = GetSquaresEventColormap(Cmap);
end
if ~exist('chansToPlot','var')
    chansToPlot = {'FZ';'CZ';'PZ';'OZ'};
end

%% Plot group results as ERPs
   
% Plot RFs
if ~isempty(group_RF)
    figure(212); clf;
    PlotResponseFnsGrid(group_RF, legendstr,tResponse,chanlocs,chansToPlot,Cmap);

    set(gcf,'Position',[1 268 568 1238]);
    for i=1:numel(chansToPlot)
        subplot(size(chansToPlot,1),size(chansToPlot,2),i);
        title('')
        ylabel(chansToPlot{i});
    end
    subplot(size(chansToPlot,1),size(chansToPlot,2),1);
    title(sprintf('%s, Group RFs',titlestr))
end

% Plot Z scores
if ~isempty(group_Z)
    figure(213); clf;
    PlotResponseFnsGrid(group_Z, legendstr,tResponse,chanlocs,chansToPlot,Cmap);

    set(gcf,'Position',[569 268 568 1238]);
    for i=1:numel(chansToPlot)
        subplot(size(chansToPlot,1),size(chansToPlot,2),i);
        title('')
        ylabel(chansToPlot{i});
    end
    subplot(size(chansToPlot,1),size(chansToPlot,2),1);
    title(sprintf('%s Group Z scores',titlestr))
end