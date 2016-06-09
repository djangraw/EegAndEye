function h = MakeLegend(marks,names,linewidths,position)

% Given the marks you'd like labeled and the names of the labels, makes a
% legend with no plot to go with it.
%
% h = MakeLegend(marks,names,linewidths,position)
%
% INPUTS:
% -marks is a cell array of length N (where N is the number of objects in 
% the plot) of strings indicating the linespec of each object in the plot,
% or an Nx3 colormap matrix (this assumes all marks are lines).
% -names is a cell array of length N of strings indicating the name of each
% object in the plot.
% -linewidths is a vector of length N indicating the linewidth of each line
% in the plot or a scalar (if they're all the same).
% -position is a 2-element vector indicating the x and y coordinate [-0,1]
% of the legend you want to create)
%
% OUTPUTS:
% -h is the handle of the axes created for the legend.
%
% Created 10/07 by DJ.
% Updated 12/18/07 by DJ.
% Updated 8/20/12 by DJ - returns axes handle.
% Updated 11/8/12 by DJ - comments
% Updated 10/2/13 by DJ - allow Nx3 colormap matrix as marks input, scalar
% as linewidths input
% Updated 9/23/14 by DJ - made it so next 'plot' command will plot to
%    original axes instead of the fake one we make in this function.


% save current axes so we can switch back to them at the end
h_original = gca;

% define position vector of axis
if nargin<3 || isempty(linewidths)
    linewidths = 1;
end
if nargin <4
    position = [0 .98]; % default position is upper right
end
position = [position .02 .02];  % plot will be small (2% of screen width/height).
if isnumeric(marks) % if it's an Nx3 colormap matrix, put each row in a cell
    marks_num = marks;
    marks = cell(1,size(marks_num,1));
    for i=1:size(marks_num,1)
        marks{i} = marks_num(i,:);
    end
end
if numel(linewidths)==1
    linewidths = repmat(linewidths,1,numel(marks));
end

% make sure number of marks are the same
if numel(marks) ~= numel(names)
    disp('ERROR: Number of marks and names must be the same.')
    return
end

% make the axes
h = axes('Position', position);
hold on

% make each mark, then put them all outside the x limits of the axis so they don't show up.
for i = 1:numel(marks)
    if isnumeric(marks{i})
        plot(1,1,'color',marks{i},'linewidth',linewidths(i))
    else
		plot(1,1,marks{i},'linewidth',linewidths(i))
    end
end
xlim([-3 0])

% make the legend, then hide the plot behind it.
legend(names)
set(gca,'visible','off')


% set current axes back to original plot (this is where the next plot
% command should plot to)
axes(h_original);
% put legend in the foreground of the figure
kids = get(gcf,'children');
set(gcf,'children',kids([2:end,1]));