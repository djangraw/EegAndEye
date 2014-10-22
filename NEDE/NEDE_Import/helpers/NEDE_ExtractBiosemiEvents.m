function all_events_cell = NEDE_ExtractBiosemiEvents(eegFilename,eventsFilename)

% Get biosemi events including port=0 and pause events.
%
% all_events_cell = NEDE_ExtractBiosemiEvents(eegFilename,eventsFilename)
%
% INPUTS:
% - eegFilename is a string indicating the .bdf file that should be read
% in using biosig's sopen function.
% - optional input eventsFilename is a string indicating where the events
% should be written to.
%
% OUTPUTS:
% - all_events_cell is an nx2 cell array, in which each row is an event.
% The first column is times (in samples) and the second is types (number
% sent to parallel port, or 'boundary' for pause events or long breaks).
%
% Created 10/20/14 by DJ.
% Updated 10/22/14 by DJ - made pause_on events into boundary events

if ~exist('eventsFilename','var')
    eventsFilename = '';
end

PAUSE_ON_CODE = 153; % determined empirically by looking at the data.
% PAUSE_OFF_CODE = 152;

% Get boundary events
dat = sopen(eegFilename);

port_events = [dat.BDF.Trigger.POS dat.BDF.Trigger.TYP-16128]; 
pauses = find(dat.BDF.Status.TYP==PAUSE_ON_CODE);
pause_events = [dat.BDF.Status.POS(pauses) dat.BDF.Status.TYP(pauses)]; 
%     resumes = find(dat.BDF.Status.TYP==PAUSE_OFF_CODE);
%     resume_events = [dat.BDF.Status.POS(resumes) dat.BDF.Status.TYP(resumes)]; 

% Add long-pause events
timebetweenevents = diff(port_events(:,1));
iStartPause = find(timebetweenevents>2500*2048/1000);
for j=1:numel(iStartPause)
    if sum(pause_events(:,1) >= port_events(iStartPause(j),1) & pause_events(:,1) <= port_events(iStartPause(j)+1,1))==0
        pause_events = [pause_events; mean(port_events(iStartPause(j)+(0:1),1)),PAUSE_ON_CODE];
    end
end
pause_events = [sort(pause_events(:,1),'ascend'),pause_events(:,2)];
sclose(dat);

% Remove any initial boundary events
if pause_events(1,1)==1
    pause_events(1,:) = [];
end

% Append all events
all_events = cat(1,port_events,pause_events); % append all events
all_events = sortrows(all_events,1); % put in chronological order

% convert pause events to boundary events and make output into cells
isPause = all_events(:,2)==PAUSE_ON_CODE;
all_events_cell = num2cell(all_events);    
all_events_cell(isPause,2) = {'boundary'};

if ~isempty(eventsFilename)
    % Write to text file
    events = cellfun(@num2str,all_events_cell,'uniformoutput',false); % convert to strings
    fileID = fopen(eventsFilename,'w');
    for i=1:size(events,1);
        fprintf(fileID,'%s %s\n',events{i,:});
    end
    fclose(fileID);
end
