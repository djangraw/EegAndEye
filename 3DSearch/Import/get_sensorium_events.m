function port_events = get_sensorium_events(filename)

% port_events = get_sensorium_events(filename)
%
% Given a sensorium filename (.dat in the current path, loads the info about
% the file and uses it to find the received parallel port events.
% - The OUTPUT port_events is an nx2 matrix, where n is the number of
% events received. The first column is the timestamp of an event (change in
% parallel port data) and the second column is the number of the event that
% was sent (see Numbers.m for a listing of event codes, which should
% correspond to the latest version of Numbers.js in the Unity project.
%
% Created 11/23/10 by DJ.

%% Get events
channels = [];
sampleRate = 1000;
nSamples = []; % autodetect length of file
offset = 0;
[~, events] = readEEG_b4preprocessing(filename,channels,nSamples,offset);


%% Convert to desired format
% arrange events into a matrix where the columns are [latency, type]
iEvents = find(diff(events)~=0)+1;
port_events = [iEvents', events(iEvents)'];