% function PlotSvdForSequenceEvents()
% Plot SVD for Squares Sequence Events
% Created 4/19/14 by DJ.

% Set up
nComponents = 3;
nEventTypes = 5;
timeWindow = [0 500];
componentSwap = [];
doPlot = 1;

%% Run
[weights, timecourse, avgW, avgTC, chanlocs, tResponse, ...
        steW, steTC,eventnames] = deal(cell(1,3));

for iExp = 1:numel(prefixes)
    prefix = prefixes{iExp};

    % Set up
    eval(sprintf('R = R_%s_sqnum;',prefix));
    % run SVD
    disp('Running SVD...')
    [weights{iExp}, timecourse{iExp}, avgW{iExp}, avgTC{iExp}, ...
        chanlocs{iExp}, tResponse{iExp}, steW{iExp}, steTC{iExp}, eventnames{iExp}] = ...
        GetGroupSvdResults(R,nComponents,nEventTypes,timeWindow,componentSwap,doPlot);
end


%% Plot
for iExp = 1:numel(prefixes)
    legendnames = eventnames{iExp};
    Cmap = GetSquaresEventColormap(legendnames); 
    figure;
    PlotGroupSvdResults(avgW{iExp}, avgTC{iExp}, ...
        chanlocs{iExp}, tResponse{iExp}, ...
        steW{iExp}, steTC{iExp}, legendnames, Cmap);
    MakeFigureTitle(sprintf('Figure %d: %s',gcf,prefixes{iExp}),0);
end

%% Stats
for iExp = 1:numel(prefixes)
    figure;
    legendnames = eventnames{iExp};
    Cmap = GetSquaresEventColormap(legendnames); 
    GetCrossSubjectStats(timecourse{iExp},tResponse{iExp},legendnames,Cmap);
    MakeFigureTitle(sprintf('Figure %d: %s',gcf,prefixes{iExp}),0);
end