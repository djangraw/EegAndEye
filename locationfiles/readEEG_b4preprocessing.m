function [eegdata,events]=readEEG_b4preprocessing(filename,channel,duration,offset, D, eventchan)

% Read in a sensorium data file and return the EEG data and events
%
% [eegdata,events]=readEEG_b4preprocessing(filename,channel,duration,offset);
%
% INPUTS:
% - filename is the name of a biosemi data file in the current path 
% (without .dat extension).
% - channel is an n-element vector of the channel numbers we want to read 
% in [1:87]
% - duration is the number of samples you want to read in starting at each 
% offset. A duration of [] indicates you want to read to the end of the
% file.
% - offset is a vector of the offsets at which you want to start reading in
% data.
%
% OUTPUTS:
% - eegdata is a matrix of the eeg data we read.  It has n rows, where n is
% the number of channels specified, and m columns, where m is the number of
% offsets specified times the duration specified.
%
% Obtained 1/8/09 from An Luo 
% Updated 11/23/10 by DJ - comments, automatic end of file checks
% Updated 2/28/13 by DJ - added D and eventchan inputs

%% Handle inputs
if nargin<2
    channel = 1:87;
end
if nargin<3
    duration = [];
end
if nargin<4
    offset = 0;
end
if nargin<5
    D=97;  % number of rows in the data file: 96 EEG channels (includes external elec's) + 1 event
end

%% Declare constants and set up
eegdata=[]; 
events=[];
blklen=2^14;  %=16384: number of samples in one 'block'

%% Open file and find end
fid=fopen([filename '.dat'],'r','b');
fseek(fid,0,'eof');
eof = ftell(fid);

%% check for errors
maxokoffset = (eof-4*4)/(4*D);
maxokduration = (eof - 4*D*max(offset) - 4*4)/(4*D);
if max(offset)>maxokoffset
    error(['You have asked the program to read beyond the file limits!\n' ...
        'Max allowable offset is %d'],maxokoffset)
elseif duration > maxokduration
    error(['You have asked the program to read beyond the file limits!\n' ...
        'given offset %d, max allowable duration is %d'],max(offset),maxokduration)
end

%% If duration was not specified, set duration to length of file
if isempty(duration)
    duration = (eof - 4*D*max(offset) - 4*4)/(4*D);
    fprintf('Setting duration to max possible given offset %d:\n duration = %d samples\n',max(offset), duration)
end

%% Read in data
for offseti=1:length(offset),
    fseek(fid,4*D*offset(offseti)+4*4,'bof');  %before  preprocessing
    %fseek(fid,4*D*offset(offseti),'bof');       %after preprocessing

    % Split the data into blocks and read it in pieces
    slices=floor(duration/blklen); % number of full blocks or 'slices' to read
    for i=1:slices
        slice = fread(fid,[D blklen],'single'); % read in 1 slice of data
        eegdata=[eegdata slice(channel+1,:)]; % append eeg data (note that first row of read-in data is events, so EEG channel i is in row i+1.)
        events=[events slice(1,:)]; % append events data
        fprintf(['Reading ' filename '.dat - offset: %d/%d - slice: %d/%d\n'],...
            offseti,length(offset),i.*blklen,duration); % display
    end
    
    % Clean up by reading the remaining data
    blklenlast=duration-(slices.*blklen); % how much data remains?
    slice =fread(fid,[D blklenlast],'single'); % read in mini-slice
    eegdata=[eegdata slice(channel+1,:)]; % append eeg data
    events=[events slice(1,:)]; % append events data
    fprintf(['Reading ' filename '.dat - offset: %d/%d - slice: %d/%d\n'],...
        offseti,length(offset),slices.*blklen+blklenlast,duration);

end

%% Clean up
fclose(fid);
