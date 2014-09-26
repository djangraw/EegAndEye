function PlotUniqueLines(x,y,marker,lineWidth,markerSize,linespergroup)

% Plots a unique set of lines so that the line styles do not repeat.
%
% PlotUniqueLines(x,y,marker,lineWidth,markerSize)
%
% INPUTS:
% - x and y are inputs to the plot command.
% - marker is a string or cell array of strings
% - lineWidth and markerSize are scalars
% - linespergroup is a scalar.
%
% Created 5/30/13 by DJ.
% Updated 8/7/13 by DJ - allow multiple markers

if isempty(x)
    x = 1:size(y,1);
end
if nargin<2
    y=x;
    x=1:size(y,1);
end
if nargin<3 || isempty(marker)
    marker = 'none';
end
if nargin<4 || isempty(lineWidth)
    lineWidth = 1;
end
if nargin<5 || isempty(markerSize)
    markerSize = 6;
end
if nargin<6 || isempty(linespergroup)
    linespergroup = size(y,2);
end

if size(x,1)==1
    x = x';
end
if size(x,2)==1
    x = repmat(x,1,size(y,2));
end

if ischar(marker)
    markers = repmat({marker},4);
else
    markers = marker;
end

% there are 7 colors, so break it down into groups of 7.
LINES_PER_GROUP = linespergroup;
Cmap = distinguishable_colors(linespergroup,{'w','k'}); % don't allow black or white
nLines = size(y,2);
nGroups = ceil(nLines/LINES_PER_GROUP);
if lineWidth==0
    lineStyles = repmat({'none'},1,4);
    lineWidth = 1;
else
    lineStyles = {'-','--',':','-.'};
end
holdState = ishold;

for i=1:nGroups
    for j=1:linespergroup
        iThis = (i-1)*LINES_PER_GROUP + j;
        if iThis>nLines
            break;
        end
        plot(x(:,iThis),y(:,iThis),'LineStyle',lineStyles{i},'Marker',markers{i},'LineWidth',lineWidth,'MarkerSize',markerSize,'Color',Cmap(j,:));
        hold on
    end
%         iThisGroup = (i-1)*LINES_PER_GROUP + (1:LINES_PER_GROUP);
%         iThisGroup(iThisGroup>nLines) = [];
%         plot(x(:,iThisGroup),y(:,iThisGroup),'LineStyle',lineStyles{i},'Marker',markers{i},'LineWidth',lineWidth,'MarkerSize',markerSize);
%         hold on
end

% Set hold state back to what it was
if holdState
    hold on
else
    hold off
end