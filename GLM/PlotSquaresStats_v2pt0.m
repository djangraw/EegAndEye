function PlotSquaresStats_v2pt0(results,eventsToPlot,iLevel)

% Plot squares GLM response functions and single-subject stats.
%
% PlotSquaresStats_v2pt0(results,eventsToPlot,iLevel)
%
% INPUTS:
% -results is a struct loaded from a RunGlmGui results struct.  It should
% have fields regressor_events, responseFns, EEG, and iLevel.
% -eventsToPlot is a vector indicating the indices of the events you want
% to plot with PlotGlmResponses.
% -iLevel is a scalar indicating the level of analysis holding the event of
% interest (default = results.iLevel).
%
% Created 7/24/12 by DJ based on PlotSquaresStats.m.
% Updated 7/26/12 by DJ - no multiple comparisons correction
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 3/19/14 by DJ - tResponse in cells

%% SACCADE-9REG (subtract out saccade regressor first)
% Handle defaults
if nargin<2 
    eventsToPlot = 1:5;
end
if nargin<3
    iLevel = [];
end

results = UpdateGlmResultsFormat(results);
if isempty(iLevel)
    iLevel = results.iLevel;
end
% Get/Plot Z scores
% [z,p] = GetGlmZscore([],results,'fdr');
% MakeFigureTitle(sprintf('%s, FDR corrected',results.filenames{results.iLevel}));
[z,p] = GetGlmZscore([],results,'none');
MakeFigureTitle(sprintf('%s, no multcompare correction',results.filenames{iLevel}));

if ~isempty(eventsToPlot)
    % Plot results
    PlotGlmResponses(results.regressor_events{iLevel}(eventsToPlot), results.responseFns{iLevel}(:,:,eventsToPlot), ...
        results.tResponse{iLevel}, results.EEG.chanlocs, [], [-6 6]);
    PlotGlmResponses(results.regressor_events{iLevel}(eventsToPlot), z(:,:,eventsToPlot), ...
        results.tResponse{iLevel}, results.EEG.chanlocs, [], [-8 8]);
    % Arrange topomovie figures
    fh=findobj(0,'type','figure');
    CascadeFigures(fh(length(eventsToPlot)*2:-1:1),length(eventsToPlot));
end