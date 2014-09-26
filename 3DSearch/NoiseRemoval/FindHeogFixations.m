function [iLeft, iRight] = FindHeogFixations(fixPos,xMin,yMax,screenCtr)

% Gets the indices of fixations to the extreme left or right of screen.
%
% [iLeft, iRight] = FindHeogFixations(fixPos,xMin,yMax,screenCtr)
%
% INPUTS:
% -fixPos is an nx2 matrix of fixation positions, where each row is the
% mean (x,y) position of one fixation.
% -xMin is a scalar indicating the minimum x position to be considered a
% fixation to the right or left.
% -yMax is a scalar indicating the maximum y position that can be 
% considered a fixation to the right or left. [default: inf (no limit)]
% -screenCenter is a 2-element vector indicating the size of the screen in
% the x and y directions. [default: [1024 768]/2]
%
% OUTPUTS:
% -iLeft is a vector of indices of fixations in fixPos that are to the left
% side of the screen.
% -iRight is a vector of indices of fixations in fixPos to the right.
%
% Created 6/4/13 by DJ.

if nargin<3 || isempty(yMax)
    yMax = inf;
end
if nargin<4 || isempty(screenCtr)
    screenCtr = [1024 768]/2;
end

doPlot = true;

% Get fixation position relative to center of screen
% fixPos = x.eyelink.fixation_positions;
fixPos = fixPos - repmat(screenCtr,size(fixPos,1),1);
% Get indices of fixations to L and R part of screen
iLeft = find(fixPos(:,1) < -xMin & abs(fixPos(:,2)) < yMax);
iRight = find(fixPos(:,1) > xMin & abs(fixPos(:,2)) < yMax);

if doPlot
    cla; hold on;
%     sacDist = x.eyelink.saccade_positions;
    plot(fixPos(:,1),fixPos(:,2),'bo');
    plot(fixPos(iLeft,1),fixPos(iLeft,2),'ro');
    plot(fixPos(iRight,1),fixPos(iRight,2),'go');
    
    xlabel('Horizontal distance')
    ylabel('Vertical distance')
    legend(sprintf('all fixations (n=%d)',length(fixPos)),...
        sprintf('leftward (n=%d)',length(iLeft)),...
        sprintf('rightward (n=%d)',length(iRight)));
end