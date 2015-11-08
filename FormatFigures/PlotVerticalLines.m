function hLines = PlotVerticalLines(values,colorspec)

% Plots vertical lines through the whole plot at the specified x locations.
%
% hLines = PlotVerticalLines(values,colorspec)
%
% INPUTS:
% -values is a vector of x positions at which to plot vertical lines.
% -colorspec is a string indicating the color and linestyle to plot.
% [default: 'k-' for a black solid line].
%
% OUTPUTS:
% -hLines is an array of handles for the lines, the same size as 'values'.
%
% Created 8/8/11 by DJ.
% Updated 4/20/15 by DJ - added hLines output.

% Handle defualts
if nargin<2 || isempty(colorspec)
    colorspec = 'k-';
end

% Plot lines
ylimits = get(gca,'YLim');
hLines = [];
for i=1:numel(values)
    if ischar(colorspec)
        hLines(i) = plot([values(i) values(i)],ylimits,colorspec);
    else
        hLines(i) = plot([values(i) values(i)],ylimits,'color',colorspec);
    end
end