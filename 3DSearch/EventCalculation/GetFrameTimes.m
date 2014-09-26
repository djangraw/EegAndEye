function frame_times = GetFrameTimes(text_file)

% Gets the times at which the camera's position was updated.
% 
% frame_times = GetFrameTimes(text_file)
%
% -input text_file should be the .asc file produced during a 3DS
% experiment.
% 
% NOTE: the Unity game settings should be set to "Sync to VBL", which makes
% each frame start just after the "vertical blank" signal of a monitor
% refresh.
%
% Created 9/20/10 by DJ.


% Setup
fid = fopen(text_file);
fseek(fid,0,'eof'); % find end of file
eof = ftell(fid);
fseek(fid,0,'bof'); % rewind to beginning

word = 'Camera at'; % the word we check for to find parallel port events
frame_times = []; %each row is timestamp, event
i = 0;          % number of messages found so far

% Get the write_ioport messages
while ftell(fid) < eof % if we haven't reached the end of the text file
    str = fgetl(fid); % read in next line of text file
    if findstr(str,word); % check for the code-word indicating a message was written
        i = i+1; % increment number of messages found
        frame_times(i,:) = sscanf(str,'MSG %d'); % read timestamp into frame_times
    end
end

% Did we find any messages?
if numel(frame_times) == 0 % if no
    warning(sprintf('Couldn''t find ''%s'' anywhere in document ''%s''!'...
        ,word,text_file))
end

% Clean up
fclose(fid);