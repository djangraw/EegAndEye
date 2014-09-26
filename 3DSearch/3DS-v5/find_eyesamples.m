function [samples, pupilsize] = find_eyesamples(text_file)

% find_eyesamples(text_file)
%
% Gets the individual eye samples from a given text file, line by line.
%
% Created 8/23/10 by DJ.

% Setup
fid = fopen(text_file);
fseek(fid,0,'eof'); % find end of file
eof = ftell(fid);
fseek(fid,0,'bof'); % rewind to beginning

word = '...'; % the word we check for to find eye samples
samples = []; % each row is time,x,y in eyelink coordinates
pupilsize = []; % size of pupil (units?)
i = 0;          % number of messages found so far

% Get the write_ioport messages
while ftell(fid) < eof % if we haven't reached the end of the text file
    str = fgetl(fid); % read in next line of text file
    if findstr(str,word); % check for the code-word indicating a message was written
        i = i+1; % increment number of messages found
        values = sscanf(str, '%*d %f %f %f'); % example line: '1234567 22.2 33.3 1111.0 ...'
        if isempty(values) % during a blink, the line will be '1234567 . . 0.0 ...'
            samples(i,1:2) = NaN;
            pupilsize(i,1) = NaN;
        else
            samples(i,1:2) = values(1:2);
            pupilsize(i,1) = values(3);
        end
    end
end

% Did we find any messages?
if numel(samples) == 0 % if no
    warning(sprintf('Couldn''t find ''%s'' anywhere in document ''%s''!'...
        ,word,text_file))
end

% Clean up
fclose(fid);