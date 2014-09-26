function [meanpos, medianpos] = EyePositionHeatMap(subject,sessions,prefix)

% [meanpos, medianpos] = EyePositionHeatMap(subject,sessions,prefix)
%
% INPUT:
% -subject is a scalar indicating the subject number.
% -sessions is an N-element vector indicating the session numbers.
% -prefix is a string indicating the experiment type code (typically '3DS', 
% 'sq', or 'sf')
%
% OUTPUT:
% -meanpos and medianpos are Nx2 matrices containing the mean and median 
% eye position (x,y) for each session.
%
% Created 2/22/13 by DJ.

%% EYE POSITION HISTOGRAM

if nargin<1
    subject = 19;
end
if nargin<2
    sessions = 2:11;
end
if nargin<3
    prefix = '3DS'; 
end


% Declare constants
n = 128; % # bins
MIN_XY = [0 0]; % (x,y) pos of bottom L of screen
MAX_XY = [1024 768]; % (x,y) pos of top R of screen
Constants = GetSquaresConstants;

% Set up
figure;
nSessions = numel(sessions); 
data = cell(1,nSessions);
nRows = ceil(sqrt(nSessions+1));
nCols = ceil((nSessions+1)/nRows);
meanpos = nan(nSessions,2);
medianpos = nan(nSessions,2);
for i=1:nSessions
    % Load
    foo = load(sprintf('%s-%d-%d-eyepos.mat',prefix,subject,sessions(i)));
    data{i} = foo.eyepos; 
    
    % Get
    [~,density,X,Y] = kde2d(data{i},n,MIN_XY,MAX_XY);

    % Plot
    subplot(nRows,nCols,i); cla; %hold on; % don't call 'hold on' until we get axis limits
    imagesc(X(1,:),Y(:,1),density);
    % superimpose stimuli
    if strcmp(prefix,'sq')
        PlotSquaresStims(Constants);
    end
    % Get mean and median eye pos
    meanpos(i,:) = nanmean(data{i},1);
    medianpos(i,:) = nanmedian(data{i},1);    
    % superimpose mean eye pos
    limits = [get(gca,'xlim'), get(gca,'ylim')]; % get axis limits
    hold on;
%     plot(meanpos(i,1),meanpos(i,2),'w+');
    plot(medianpos(i,1),medianpos(i,2),'r+','markersize',10);
    axis(limits); % put limits back
    axis equal
    % annotate plot
    colorbar
    set(gca,'clim',[0 2e-4],'ydir','reverse');
    title(sprintf('%s-%d-%d',prefix,subject,sessions(i)));    
end

% Average
subplot(nRows,nCols,numel(sessions)+1); cla;
[bandwidth,density,X,Y] = kde2d(cat(1,data{:}),n,MIN_XY,MAX_XY);
imagesc(X(1,:),Y(:,1),density);
if strcmp(prefix,'sq')
    PlotSquaresStims(Constants);
end
colorbar
set(gca,'clim',[0 2e-4]);
title(sprintf('%s-%d, all %d sessions',prefix,subject,numel(sessions)));

end % function TEMP_MapEogVoltages


function PlotSquaresStims(Constants)
    hold on
    plot(Constants.LEFTCROSS_X, Constants.LEFTCROSS_Y,'w+','markersize',10);
    plot(Constants.SQUARE_X, Constants.SQUARE_Y,'ws','markersize',10);
    plot(Constants.RIGHTCROSS_X, Constants.RIGHTCROSS_Y,'w+','markersize',10);    
end