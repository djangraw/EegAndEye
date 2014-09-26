function ImageLevel(levelname, offset)

% Plot an image with a given offset.
%
% ImagePoints(Cimage, offset)
%
% INPUTS:
% -levelname is a string that indicates an image filename in the current
% path.
% - offset is the x,y position of the origin in the image.  It will be
% subtracted from the x and y axes when plotting the image.
%
% Created 4/18/11 by DJ.

% handle inputs
if nargin<2
    offset = [0 0];
end

C = imread(levelname);

% display image
imagesc((1:size(C,2))-offset(1), (1:size(C,1))-offset(2), C);
% format
axis image;
set(gca,'ydir','normal');
hold on;
% make title
title(levelname);
