function [timestamps endpos] = find_saccadeevents(text_file)

% [timestamps endpos] = find_saccadeevents(text_file);
%
% - This function takes an eyelink text file and extracts the times when
% saccades took place.
% - The INPUT text_file should be the filename of the eyelink text file.
% - To create the eyelink text file, run the .edf file created during a
% 3D Search (3DS) experiment through the program Visual EDF2ASC.
% - The OUTPUT timestamps will be an n-element vector, where n is the 
% number of saccades.  Each element is an EyeLink timestamp of the end time
% of a saccade.
% - The OUTPUT endpos will be an nx2 matrix, where n is the number of
% saccades.  Each row is the x and y position of the end of the saccade.
% NOTE: don't forget to recalibrate using the values found in Unity.
%
% Created 6/14/10 by DJ.
% Updated 6/17/10 by DJ - comments
% Updated 6/24/10 by DJ - use Saccade reports instead of messages
% Updated 8/25/10 by DJ - made timestamps a column vector
% Updated 9/13/10 by DJ - comments
% Updated 10/15/10 by DJ - added check for proper saccade logging

% Setup
fid = fopen(text_file);
fseek(fid,0,'eof'); % find end of file
eof = ftell(fid);
fseek(fid,0,'bof'); % rewind to beginning

word = 'ESACC'; % the word we check for to find saccade events
timestamps = []; % start with an empty array of timestamps
i = 0;          % number of elements currently in timestamps


% Get the Timestamps
while ftell(fid) < eof % if we haven't found our saccade marker
    str = fgetl(fid); % read in another line from the text file
    if findstr(str,word); % if this line contains a saccade timestamp
        values = sscanf(str,'ESACC %*s %*d %d %*d %*f %*f %f %f'); % extract the saccade end time and end position from this regularly formatted line:
        if numel(values)==3
            i = i+1; % increment the number of timestamps found
            % (ESACC <eye> <starttime> <endtime> <duration> <startx> <starty> <endx> <endy> <amp> <peakvelocity>)        
            timestamps(i,1) = values(1);
            endpos(i,1:2) = values(2:3);
        else
            warning('Saccade End Logged Improperly: %s', str);
        end
    end
end

% Did we find the word anywhere?
if numel(timestamps) == 0 % if no
    warning(sprintf('Couldn''t find ''%s'' anywhere in document ''%s''!'...
        ,word,text_file))
end

% Clean up
fclose(fid);