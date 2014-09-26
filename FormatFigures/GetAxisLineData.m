function [data,colors,displaynames] = GetAxisLineData(hAxis)

% Extracts data & info about the line plots on the current axis.
%
% [data,colors,displaynames] = GetAxisLineData(hAxis)
%
% INPUTS:
% -hAxis is a handle to an axis [default: gca]
%
% OUTPUTS:
% -data is an N-element vector of cells, where N is the number of line
% plots on the current axis. Each cell contains a 2xM or 3xM matrix in 
% which the first row is the x data, the 2nd row is the y data, and the 3rd 
% row is the z data (if it exists).
% -colors is an N-element vector of cells containing the colors of each
% line (generally an RGB vector).
% -displaynames is an N-element vector of cells containing the names of
% each line (if the axis has a legend).
%
% SAMPLE USAGE:
% To make a random line plot and then recreate it:
% % Make a random plot in a new figure
% figure;
% plot(rand(5));
% legend('one','two','three','four','five');
% % Extract the data
% [data, colors, displaynames] = GetAxisLineData(gca);
% % Remake the plot in another figure
% figure; cla; hold on;
% for i=1:numel(data)
%     plot(data{i}(1,:),data{i}(2,:),'Color',colors{i})
% end
% legend(displaynames)
%
% Created 5/8/13 by DJ.

if nargin<1
    hAxis = gca;
end

% Set up
kids = get(hAxis,'Children');
types = get(kids,'Type');
isLine = strcmp('line',types);
linekids = flipud(kids(isLine)); % flip to put lines in legend order

% Get colors and names
colors = get(linekids,'Color');
displaynames = get(linekids,'DisplayName');
% Get data
xdata = get(linekids,'XData');
ydata = get(linekids,'YData');
zdata = get(linekids,'ZData');

% Combine data
N = numel(linekids);
data = cell(1,N);
for i=1:N
   if ~isempty(zdata{i})
       data{i} = [xdata{i}; ydata{i}; zdata{i}];
   else
       data{i} = [xdata{i}; ydata{i}];
   end
end