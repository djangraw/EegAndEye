function EEG = NEDE_AddHeogEvents(EEG0, y, xMin, yMax)

% Adds fixation start/end events for saccades to the sides of the screen.
% 
% EEG = NEDE_AddHeogEvents(EEG0, y, xMin, yMax)
%
% INPUTS:
% -EEG0 is an eeglab struct made from a combined n sessions.
% -y is an n-element vector of NEDE structs for the corresponding n 
% sessions, each including an 'events' field with fixation information.
% -xMin is a scalar indicating the minimum x position to be considered a
% fixation to the right or left. [default: 200]
% -yMax is a scalar indicating the maximum y position that can be 
% considered a fixation to the right or left. [default: 50]
%
% OUTPUTS:
% -EEG is the same as EEG0, but with 'FSL' & 'FEL' events for fixation 
% starts & ends on the left side of the screen, and 'FSR' & FER' for those
% on the right side.
%
% Created 9/26/14 by DJ based on AddHeogFixationMarkers.m.


% Set defaults
if nargin<3 || isempty(xMin)
    xMin = 200;
end
if nargin<4 || isempty(yMax)
    yMax = 50;
end
screenCenter = [y(1).params.screen.width, y(1).params.screen.height]/2;

% Declare variables
nSessions = numel(y);
allEvents = cell(1,nSessions);
allCodes = cell(1,nSessions);
figure;
nRows = ceil(sqrt(nSessions));
nCols = ceil(nSessions/nRows);

% Get HEOG saccades
for i=1:nSessions
    % Set up plot made in FindHeogFixations
    subplot(nRows,nCols,i);

    % Pick out fixations    
    [iLeft, iRight] = FindHeogFixations(y(i).events.fixation.position,xMin,yMax,screenCenter);

    % Get saccade times       
    startTimesLeft = y(i).events.fixation.time_start(iLeft);
    endTimesLeft = y(i).events.fixation.time_end(iLeft);
    startTimesRight = y(i).events.fixation.time_start(iRight);
    endTimesRight = y(i).events.fixation.time_end(iRight);

        
    % Make events matrix
    allEvents{i} = cat(1,startTimesLeft,endTimesLeft,startTimesRight,endTimesRight);
    allCodes{i} = cat(1,repmat({'FSL'},size(startTimesLeft)),...
        repmat({'FEL'},size(endTimesLeft)), ...
        repmat({'FSR'},size(startTimesRight)), ...
        repmat({'FER'},size(endTimesRight)));
    
end

% Add events to EEG struct
EEG = NEDE_AddEeglabEvents(y,EEG0,allEvents,allCodes);