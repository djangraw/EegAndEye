function ImagePoints(points, offset, color, marker, markersize)

% Plot an image with points overlaid on it.
%
% ImagePoints(points, offset, color, marker)
%
% INPUTS:
% - points is a kx2 matrix where each row is the x,y position of a point in
% pixel coordinates.
% - offset is the x,y position of the origin in the image.  It will be
% subtracted from the x and y axes when plotting the image. [default: 0 0]
% - color is a colorspec indicator (either a string like 'g' or a 3-element
% rgb vector like [0 0.4 0.6]).
% - marker is a plot/scatter mark string like 'o' or '.'
% -markersize is a scalar value indicating how big the markers should be.
%
% Created 4/18/11 by DJ.

% handle inputs
if nargin<2
    offset = [0 0];
end
if nargin<3
    color = 'g';
end
if nargin<4
    marker = '.';
end
if nargin<5
    markersize = 10;
end

% offset points
points_offset = points + repmat(offset,size(points,1),1);    

% plot points
if ~isempty(points_offset)
    plot(points_offset(:,1),points_offset(:,2),marker,'MarkerEdgeColor',color,'MarkerSize',markersize);
end