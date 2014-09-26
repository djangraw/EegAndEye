function PlotResponseFnsGrid(responseFns,regressor_events,tResponse,chanlocs,gridOfChans,colors)

% PlotResponseFnsGrid(responseFns,regressor_events,tResponse,chanlocs,gridOfChans,colors)
%
% Created 3/8/13 by DJ.
% Updated 10/8/13 by DJ - added colors input

if ~exist('colors','var')
    colors = [];
end

% Set up
[nRows, nCols] = size(gridOfChans);
clf;
isEmptyPlot = false(nRows,nCols);
% Plot
for i=1:nRows
    for j=1:nCols
        % Make plot
        iPlot = (i-1)*nCols + j;
        subplot(nRows,nCols,iPlot);
        % Find channel or electrode name to plot
        if iscell(gridOfChans)
            chanToPlot = gridOfChans{i,j};
        else
            chanToPlot = gridOfChans(i,j);
        end
        if isempty(chanToPlot) % to decide where to put legends later
            isEmptyPlot(i,j) = true;
        end
        % Plot response functions
        PlotResponseFns(responseFns,regressor_events,tResponse,chanlocs,chanToPlot,colors)
    end
end

% Standardize y scale
max_ylim = [0 0];
for iPlot = 1:(nRows*nCols)
    this_ylim = get(subplot(nRows,nCols,iPlot),'ylim');
    max_ylim(1) = min(max_ylim(1),this_ylim(1));
    max_ylim(2) = max(max_ylim(2),this_ylim(2));
end
for iPlot = 1:(nRows*nCols)    
    set(subplot(nRows,nCols,iPlot),'ylim',max_ylim);   
    plot([0 0],get(gca,'ylim'),'k')
end

% Make sure there's only one legend
if ~any(isEmptyPlot(:))
    isEmptyPlot(1,1) = true;
end
for i=1:nRows
    for j=1:nCols
        % Make plot
        iPlot = (i-1)*nCols + j;
        if ~isEmptyPlot(i,j)
            legend(subplot(nRows,nCols,iPlot),'hide');    
        end
    end
end
