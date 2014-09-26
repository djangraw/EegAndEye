function PlotVerticalLines(values,colorspec)

% Plots vertical lines through the whole plot at the specified x locations.
%
% PlotVerticalLines(values,colorspec)
%
% INPUTS:
% -values is a vector of x positions at which to plot vertical lines.
% -colorspec is a string indicating the color and linestyle to plot.
% [default: 'k-' for a black solid line].
%
% Created 8/8/11 by DJ.

% Handle defualts
if nargin<2 || isempty(colorspec)
    colorspec = 'k-';
end

% Plot lines
ylimits = get(gca,'YLim');
for i=1:numel(values)
    plot([values(i) values(i)],ylimits,colorspec);
end