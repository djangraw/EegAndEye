function hNew = CopyAxisToNewFigure(hAxis,fignum)

% Created 10/1/13 by DJ.

% Make new figure
figure(fignum);
% Get standard axis position
clf; axes;
stdpos = get(gca,'Position');
clf;

hNew = copyobj(hAxis,fignum);
set(hNew,'Position',stdpos);
