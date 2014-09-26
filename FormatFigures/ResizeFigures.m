function ResizeFigures(figs,newSize,sizeofmatrix)

% Make all figures the specified size.
%
% ResizeFigures(figs,newSize,sizeofmatrix)
%
% Created 3/13/14 by DJ.

if ~exist('newSize','var') || isempty(newSize)
    newSize = [inf inf];
end

if exist('sizeofmatrix','var')
    cascade = true;
else
    cascade = false;
end

nFigs = numel(figs);
% Get size
[figPos,figSize] = deal(zeros(nFigs,2));
for i=1:nFigs
    foo = get(figs(i),'position');
    figPos(i,:) = foo(1:2);
    figSize(i,:) = foo(3:4);
end


% Get width and height of screen in pixels
% NOTE: menuBarHeight (pixels) is the height of the menu bar at the top of
% each window.  It's determined empirically by dragging a figure to the 
% top of the screen, getting its y position, and subtracting it from the 
% screen height. If it's too small, the program will be slow.  If it's too
% large, the first row of figures won't be at the top of the screen.
% DockHeight (pixels) is the size of the dock at the bottom of the screen.
% It's determined empirically by manually resizing a window to take up the
% whole screen (vertically), getting its height, and subtracting that from
% the screen height minus the menu bar height.
menuBarHeight = 100;
dockHeight = 100;
set(0,'Units','pixels');
screenSize = get(0,'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4)-menuBarHeight-dockHeight;

% Adjust size to fit within maximum allowed
if cascade
    maxSize = [screenWidth/sizeofmatrix(2), screenHeight/sizeofmatrix(1)];
else
    maxSize = screenSize(3:4);
end
newSize = min(newSize,maxSize);

% Resize
for i=1:nFigs
    set(figs(i),'position',[figPos(i,:),newSize]);
end

% Cascade
if cascade
    CascadeFigures(figs,sizeofmatrix(2));
end

