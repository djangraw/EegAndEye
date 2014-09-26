function [ts_events] = find_writeioport(text_file)

% [ts_events] = find_writeioport(text_file);
%
% - This function takes an eyelink text file and extracts the times when
% events were sent via parallel port to the EEG.  This can be compared
% against the times when the EEG received these events.
% - The INPUT text_file should be the filename of the eyelink text file.
% - To create the eyelink text file, run the .edf file created during a
% 3D Search (3DS) experiment through the program Visual EDF2ASC.
% - The OUTPUT ts_events will be an nx2 matrix, where n is the number of
% events.  The first column is the timestamp of the event, and the second
% column is the event number (see Numbers.js for the meaning of these
% codes).
%
% Created 6/14/10 by DJ.
% Updated 6/17/10 by DJ - comments
% Updated 8/18/10 by DJ - eliminated starttime.

% Setup
fid = fopen(text_file);
fseek(fid,0,'eof'); % find end of file
eof = ftell(fid);
fseek(fid,0,'bof'); % rewind to beginning

word = 'write_ioport'; % the word we check for to find parallel port events
ts_events = []; %each row is timestamp, event
i = 0;          % number of messages found so far

% Get the write_ioport messages
while ftell(fid) < eof % if we haven't reached the end of the text file
    str = fgetl(fid); % read in next line of text file
    if findstr(str,word); % check for the code-word indicating a message was written
        i = i+1; % increment number of messages found
        ts_events(i,:) = sscanf(str,'MSG %d !CMD %*d write_ioport 0x378 %d'); % read timestamp and message numbers into ts_events
    end
end

% Did we find any messages?
if numel(ts_events) == 0 % if no
    warning(sprintf('Couldn''t find ''%s'' anywhere in document ''%s''!'...
        ,word,text_file))
end

% Clean up
fclose(fid);