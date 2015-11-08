function hLines = PlotHorizontalLines(values,colorspec)

% Plots horizontal lines through the whole plot at the specified y locations.
%
% hLines = PlotHorizontalLines(values,colorspec)
%
% INPUTS:
% -values is a vector of x positions at which to plot vertical lines.
% -colorspec is a string indicating the color and linestyle to plot.
% [default: 'k-' for a black solid line].
%
% OUTPUTS:
% -hLines is an array of handles for the lines, the same size as 'values'.
%
% Created 10/20/14 by DJ based on PlotVerticalLines.
% Updated 4/20/15 by DJ - added hLines output.

% Handle defualts
if nargin<2 || isempty(colorspec)
    colorspec = 'k-';
end

% Plot lines
xlimits = get(gca,'XLim');
hLines = [];
for i=1:numel(values)
    hLines(i) = plot(xlimits,[values(i) values(i)],colorspec);
end