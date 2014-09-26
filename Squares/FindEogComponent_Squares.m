function [component, offset] = FindEogComponent_Squares(EEG,y,elecToPlot)

% Created 6/28/12 by DJ.

if nargin<3
    elecToPlot = 'F8';
end

% Get pixel coordinates of center of screen
Constants = GetSquaresConstants;
X_CENTER = Constants.SQUARE_X(3);
Y_CENTER = Constants.SQUARE_Y(3);

nSessions = numel(y);
boundaryTimes = [0 EEG.event(strmatch('boundary', {EEG.event(:).type})).latency]; % in seconds
if numel(boundaryTimes)<nSessions-1
    error('%d sessions given as input but only %d boundary events!',...
        nSessions,numel(boundaryTimes))
end

fixStartTime = cell(1,nSessions);
fixEndTime = cell(1,nSessions);

% get relevant info
for i=1:nSessions
    % Set up syncs
    elSyncs = y(i).sync.eyelink;
    eegSyncs = y(i).sync.eeg + boundaryTimes(i); % convert to EEG samples
    % Convert eyelink to eeg times
    fixStartTime{i} = EyelinkToEegTimes(y(i).fixation.start_time,elSyncs,eegSyncs);
    fixEndTime{i} = EyelinkToEegTimes(y(i).fixation.end_time,elSyncs,eegSyncs);
    % Get position info
    xpos{i} = y(i).fixation.position(:,1) - X_CENTER;
    ypos{i} = y(i).fixation.position(:,2) - Y_CENTER;
    sqnum{i} = y(i).fixation.squarenum;
end
% convert to big long vectors
fixStartAll = cat(1,fixStartTime{:});
fixEndAll = cat(1,fixEndTime{:});
xposAll = cat(1,xpos{:});
yposAll = cat(1,ypos{:});
sqnumAll = cat(1,sqnum{:});
nFixations = length(fixStartAll);


% Get mean data in every fixation
meanData = nan(EEG.nbchan,nFixations);
for i=1:nFixations
    EEGsamples = ceil(fixStartAll(i)):floor(fixEndAll(i));
    if all(EEGsamples<size(EEG.data,2))
        meanData(:,i) = mean(EEG.data(:,EEGsamples),2);
    end
end

% Group fixations to each square
groupMedian = nan(EEG.nbchan,7);
groupXpos = nan(1,7);
for j=1:7
    isInGroup = sqnumAll==j-1;
    groupMedian(:,j) = median(meanData(:,isInGroup),2);
    groupXpos(j) = median(xposAll(isInGroup));
end

% Get horizontal EOG component
isInSquare = ~isnan(sqnumAll);
for i=1:EEG.nbchan
    b = regress(groupMedian(i,:)', [groupXpos', ones(7,1)]);
    component(i) = b(1);
    offset(i) = b(2);
end
% betas = regress(groupXpos',[groupMedian', ones(7,1)]);
% component = 1./betas(1:end-1);
% offset = -betas(end)./betas(1:end-1);

% Plot regression line
figure(555)
iElec = find(strcmp(elecToPlot,{EEG.chanlocs.labels}));
subplot(2,4,1);
cla; hold on;
for j=1:7
    isInGroup = sqnumAll==j-1;
    scatter(xposAll(isInGroup),meanData(iElec,isInGroup),'o');
end
scatter(groupXpos,groupMedian(iElec,:),'k+');
xlim = [min(xposAll(~isnan(sqnumAll))), max(xposAll(~isnan(sqnumAll)))];
plot(xlim, xlim*component(iElec)+offset(iElec))
ylabel('mean value on this electrode (uV)');
xlabel('x position of eye');
title(sprintf('EOG regression for electrode %s',EEG.chanlocs(iElec).labels))
% legend('mean EEG during fixation', 'fit');

% Plot channel values
for j=1:7
    subplot(2,4,j+1);
    topoplot(groupMedian(:,j),EEG.chanlocs,'electrodes','on');
    colorbar;
    title(sprintf('Position %d',j-1));
end

% Plot component & offset
figure(556)
subplot(1,2,1)
topoplot(component,EEG.chanlocs,'electrodes','on');
title('HEOG component')
colorbar
subplot(1,2,2)
topoplot(offset,EEG.chanlocs,'electrodes','on');
title('HEOG offset')
colorbar



