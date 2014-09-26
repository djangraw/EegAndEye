function [eegdata, events] = readEEG_sensorium2013(filename,eegchan,D,eventchan)

% Read in a sensorium data file and return the EEG data and events
%
% [eegdata,events]=readEEG_raw(filename,eegchan,D,eventchan);
%
% INPUTS:
% - filename is the name of a biosemi data file in the current path 
% (without .dat extension).
% - eegchan is an n-element vector of the channel numbers we want to read 
% in [9:92]
% - D is a scalar indicating the number of channels (found in the header
% file).
% - eventchan is a scalar indicating the channel numbers of the events.
%
% OUTPUTS:
% - eegdata is a matrix of the eeg data we read.  It has n rows, where n is
% the number of channels specified, and m columns, where m is the number of
% offsets specified times the duration specified.
%
% Created 2/28/13 by DJ based on readEEG_b4preprocessing.m
% Updated 4/23/13 by DJ - defaults

%% Handle inputs
if nargin<2 || isempty(eegchan)
    eegchan = 9:92; % exclude 2 digital channels, 6 eog channels
end
if nargin<3 || isempty(D)
    D = 98; % number of rows in the data file: 96 EEG channels (includes external elec's) + 1 event
end
if nargin<4 || isempty(eventchan)
    eventchan = 1;
end

%% Read in data
fprintf(['Reading ' filename '.dat...\n']); % display
fid=fopen([filename '.dat'],'r','b'); % Open file
alldata = fread(fid,[D inf],'single'); % read in all the data
eegdata=alldata(eegchan,:); % append eeg data (note that first row of read-in data is events, so EEG channel i is in row i+1.)
events=alldata(eventchan,:); % append events data

%% Clean up
fclose(fid);

fprintf('Done!\n')