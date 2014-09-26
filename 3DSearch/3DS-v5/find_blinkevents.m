function blink_times = find_blinkevents(text_file)

% blink_times = find_blinkevents(text_file)
%
% - This function takes an eyelink text file and extracts the times when
% the subject blinked.  
% - The INPUT text_file should be the filename of the eyelink text file.
% - To create the eyelink text file, run the .edf file created during a
% 3D Search (3DS) experiment through the program Visual EDF2ASC.
% - The OUTPUT button_times will be an n-element vector, where each element is
% the eyelink time at which the eyelink system detected a blink.
%
% Created 10/18/10 by DJ (based on find_writeioport)


% Setup
fid = fopen(text_file);
fseek(fid,0,'eof'); % find end of file
eof = ftell(fid);
fseek(fid,0,'bof'); % rewind to beginning

word = 'EBLINK'; % the word we check for to find relevant events
blink_times = []; %each row is timestamp, event
i = 0;          % number of messages found so far

% Get the write_ioport messages
while ftell(fid) < eof % if we haven't reached the end of the text file
    str = fgetl(fid); % read in next line of text file
    if findstr(str,word) % check for the code-word indicating a message was written
        values = sscanf(str,'EBLINK %*s %d %*d'); % Message format: EBLINK <eye (R/L)> <blinkstart> <blinkend>
        if numel(values)>0
            i = i+1; % increment number of messages found
            blink_times(i) = values(1); % read timestamp and message numbers into ts_events
        end
    end
end

% Did we find any messages?
if i == 0 % if no
    warning(sprintf('Couldn''t find ''%s'' anywhere in document ''%s''!'...
        ,word,text_file))
end

% Clean up
fclose(fid);