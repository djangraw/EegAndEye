function CascadeFigures(figNumbers,nCols)

% Takes the specified figures and places them in orderly rows and columns.
%
% CascadeFigures(figNumbers,nCols)
%
% INPUTS:
% -figNumbers is a vector of integers indicating the numbers of the figures
% you want to reposition.
% -nCols is a scalar indicating the number of columns you want the figures
% to show up in (that is, the number of figures per row). If nCols=0, this
% indicates a special case where you want all the figures in the upper-left
% corner.
%
% Note: parameter menuBarHeight may be different for each machine.  The
% wrong value could make this program run slow or position the figures
% below the top of the screen.  Edit this value as specified in the code.
%
% Created 6/27/12 by DJ.
% Updated 11/13/12 by DJ - added nCols=0 special case

% Get width and height of screen in pixels
% NOTE: menuBarHeight (pixels) is the height of the menu bar at the top of
% each window.  It's determined empirically by dragging a figure to the 
% top of the screen, getting its y position, and subtracting it from the 
% screen height. If it's too small, the program will be slow.  If it's too
% large, the first row of figures won't be at the top of the screen.
menuBarHeight = 21; 
set(0,'Units','pixels');
screenSize = get(0,'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4)-menuBarHeight;

% Special case: if nCols = 0, put all figures in top left
if nCols==0
    oldPos = get(figure(figNumbers(1)),'OuterPosition');
    newPos = [0 screenHeight oldPos([3 4])];
    for iFig=1:numel(figNumbers)
        set(figure(figNumbers(iFig)),'OuterPosition',newPos);
    end
    return;
end

% Calculate number of rows & columns in figure array
nFigs = numel(figNumbers);
nRows = ceil(nFigs/nCols);

% Get width and height of each figure
width = zeros(1,nRows*nCols);
height = zeros(1,nRows*nCols);
for iFig=1:nFigs
    pos = get(figure(figNumbers(iFig)),'OuterPosition');
    width(iFig) = pos(3);
    height(iFig) = pos(4);
end

% Calcluate row heights
rowHeight = zeros(1,nRows);
for i=1:nRows
    inRow = (i-1)*nCols + (1:nCols);
    rowHeight(i) = max(height(inRow));
end
if sum(rowHeight)>screenHeight
    rowHeight(:) = screenHeight/nRows;
end

% Calculate column widths
colWidth = zeros(1,nRows);
for j=1:nCols
    inCol = (0:nCols:nFigs-1) + j;
    colWidth(j) = max(width(inCol));
end
if sum(colWidth)>screenWidth
    colWidth(:) = screenWidth/nCols;
end
        
% Position each figure
for i=1:nRows           
    for j=1:nCols
        % Calculate figure number
        iFig = (i-1)*nCols + j;
        if iFig>nFigs 
            break; 
        end
        % Get new position
        newPos = [sum(colWidth(1:j-1)), screenHeight-sum(rowHeight(1:i)), width(iFig), height(iFig)];
        % Move figure
        set(figure(figNumbers(iFig)),'OuterPosition',newPos);
    end
end
        