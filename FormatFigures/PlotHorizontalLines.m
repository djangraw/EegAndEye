function PlotHorizontalLines(values,colorspec)

% Plots horizontal lines through the whole plot at the specified y locations.
%
% PlotHorizontalLines(values,colorspec)
%
% INPUTS:
% -values is a vector of x positions at which to plot vertical lines.
% -colorspec is a string indicating the color and linestyle to plot.
% [default: 'k-' for a black solid line].
%
% Created 10/20/14 by DJ based on PlotVerticalLines.

% Handle defualts
if nargin<2 || isempty(colorspec)
    colorspec = 'k-';
end

% Plot lines
xlimits = get(gca,'XLim');
for i=1:numel(values)
    plot(xlimits,[values(i) values(i)],colorspec);
end