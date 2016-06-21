function [xData, yData] = GetBarPositions(hBar)

% xData = GetBarPositions(hBar)
%
% For each set of bars in a barplot, find and return the centers of the
% bars.
%
% INPUTS:
% -hBar is an N-element vector of handles to a barplot.
% 
% OUTPUTS:
% -xData is an NxM matrix of the x positions of each bar.
%
% Created 6/21/16 by DJ.

drawnow; %allows the figure to be created
[xData,yData] = deal(nan(numel(hBar),numel(hBar(1).XData)));
for iBar = 1:numel(hBar)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData(iBar,:) = hBar(iBar).XData+hBar(iBar).XOffset;
    yData(iBar,:) = hBar(iBar).YData;
end