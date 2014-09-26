function [objects, objectnames, objectlocations, objecttimes, objectisessions] = GetObjectList(subject,sessions,eyelinkoreeg)

% Make a list of the names of all the objects used in all sessions.
%
% [objects, objectnames, objectlocations] = GetObjectList(subject,sessions)
% 
% INPUTS:
% -subject is the subject number and sessions is a vector of sessions
% numbers of 3DSearch data files as imported with Import_3DS_Data_v3.
% -eyelinkoreeg is a string ['eyelink' or 'eeg'] indicating which device
% you would like to provide the times for objecttimes.
%
% OUTPUTS:
% -objects is an array of structs that is all the x.object structs
% concatenated.
% -objectnames is a cell array of strings indicating the object name of each
% object placed in each session.
% -objectlocations is an array of the location of each object.
%
% Created 4/18/11 by DJ.
% Updated 11/7/11 by DJ - added eyelinkoreeg option, objectsessions output.
% Updated 12/1/11 by DJ - changed objectisessions output (session INDEX)

% set defaults
if nargin<3
    eyelinkoreeg = 'eeg';
end

% initialize
objects = [];
objectnames = [];
objectlocations = [];
objecttimes = [];
objectisessions = [];
GetNumbers;
% main loop
for i = 1:numel(sessions)
    % load
    load(sprintf('3DS-%d-%d',subject,sessions(i)));
    eventTimes = x.(eyelinkoreeg).object_events(:,1); % convert to units of seconds
    eventObjects = x.(eyelinkoreeg).object_events(:,2)-Numbers.ENTERS;
    % add to arrays
    objects = [objects x.objects]; % all object info
    objectnames = [objectnames {x.objects.name}]; % object names ('<category>-<number>')
    objectisessions = [objectisessions repmat(i,1,numel(x.objects))]; % INDEX of session, not session number
    for j=1:numel(x.objects) 
        % get x,z location
        objectlocations = [objectlocations; ...
            x.objects(j).createposition(1), x.objects(j).createposition(3)];
        appearEvent = find(eventObjects == j);
        % get appear time if it appeared
        if isempty(appearEvent)
            objecttimes = [objecttimes, NaN];
        else
            objecttimes = [objecttimes, eventTimes(appearEvent)];
        end
    end
    
    
end
switch eyelinkoreeg
    case 'eyelink'
        objecttimes = objecttimes/1000;
    case 'eeg'
        objecttimes = objecttimes/x.eeg.eventsamplerate;
end