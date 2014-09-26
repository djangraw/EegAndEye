function EEG = AddHeogFixationMarkers(EEG0, y, xMin, yMax, screenCenter)

% Adds fixation start/end events for saccades to the sides of the screen.
% 
% EEG = AddHeogFixationMarkers(EEG0, y, xMin, yMax, screenCenter)
%
% INPUTS:
% -EEG0 is an eeglab struct made from a combined n sessions.
% -y is an n-element vector of behavior structs for the corresponding n 
% sessions, each including an 'eyelink' field with fixation information.
% -xMin is a scalar indicating the minimum x position to be considered a
% fixation to the right or left. [default: 200]
% -yMax is a scalar indicating the maximum y position that can be 
% considered a fixation to the right or left. [default: 50]
% -screenCenter is a 2-element vector indicating the size of the screen in
% the x and y directions. [default: [1024 768]/2]
%
% OUTPUTS:
% -EEG is the same as EEG0, but with 'FSL' & 'FEL' events for fixation 
% starts & ends on the left side of the screen, and 'FSR' & FER' for those
% on the right side.
%
% Created 6/4/13 by DJ.
% Updated 3/10/14 by DJ - adapted to new squares data format

% Set defaults
if nargin<3 || isempty(xMin)
    xMin = 200;
end
if nargin<4 || isempty(yMax)
    yMax = 50;
end
if nargin<5 || isempty(screenCenter)
    screenCenter = [1024 768]/2;
end

% % Load data
% fprintf('loading data...\n');
% eegFilename = sprintf('%s-%d-all%s.set',prefix,subject,eegSuffix);
% EEG = pop_loadset(eegFilename);
% y = loadBehaviorData(subject,sessions,prefix);

% Declare variables
nSessions = numel(y);
allEvents = cell(1,nSessions);
allCodes = cell(1,nSessions);
figure;
nRows = ceil(sqrt(nSessions));
nCols = ceil(nSessions/nRows);

% Get HEOG saccades
for i=1:nSessions
    % Get and plot saccades
    subplot(nRows,nCols,i);
    if isfield(y(i),'eyelink')
        
        [iLeft, iRight] = FindHeogFixations(y(i).eyelink.fixation_positions,xMin,yMax,screenCenter);

        % Get saccade times       
        startTimesLeft = y(i).eyelink.fixation_times(iLeft,1);
        endTimesLeft = y(i).eyelink.fixation_times(iLeft,2);
        startTimesRight = y(i).eyelink.fixation_times(iRight,1);
        endTimesRight = y(i).eyelink.fixation_times(iRight,2);
    else
        [iLeft, iRight] = FindHeogFixations(y(i).fixation.position,xMin,yMax,screenCenter);

        % Get saccade times       
        startTimesLeft = y(i).fixation.start_time(iLeft);
        endTimesLeft = y(i).fixation.end_time(iLeft);
        startTimesRight = y(i).fixation.start_time(iRight);
        endTimesRight = y(i).fixation.end_time(iRight);
    end
        
    % Make events matrix
    allEvents{i} = cat(1,startTimesLeft,endTimesLeft,startTimesRight,endTimesRight);
    allCodes{i} = cat(1,repmat({'FSL'},size(startTimesLeft)),...
        repmat({'FEL'},size(endTimesLeft)), ...
        repmat({'FSR'},size(startTimesRight)), ...
        repmat({'FER'},size(endTimesRight)));
    
end

% Add events to EEG struct
y = AddSyncField(y);
EEG = AddEeglabEvents_MultiSession(EEG0,y,allEvents,allCodes);