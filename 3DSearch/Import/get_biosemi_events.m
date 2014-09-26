function port_events = get_biosemi_events(filename)

% port_events = get_biosemi_events(filename)
%
% Given a biosemi filename (.bdf) in the current path, loads the info about
% the file and uses it to find the received parallel port events.
% - The OUTPUT port_events is an nx2 matrix, where n is the number of
% events received. The first column is the timestamp of an event (change in
% parallel port data) and the second column is the number of the event that
% was sent (see Numbers.m for a listing of event codes, which should
% correspond to the latest version of Numbers.js in the Unity project.
%
% Created 6/16/10 by DJ.

% Open file
dat = sopen(filename);
%DAT = sread(dat, Inf);
%EEG = biosig2eeglab(dat,DAT,[],[]);

% Get event sample number (POS) and number (TYP)
% For some reason, all events are offset by 16128 (dec2hex('3F00'))
port_events = [dat.BDF.Trigger.POS dat.BDF.Trigger.TYP-16128]; 

% Clean up
sclose(dat);