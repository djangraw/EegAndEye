function varargout = find_events_lb(text_file,event_type)

% [varargout] = find_events(text_file,event_type)
%
% EXAMPLE USAGE:
% blink_start_times = find_events(text_file,'blink');
% button_press_times = find_events(text_file,'button');
% [eye_xy, pupil_size] = find_events(text_file,'eyesample');
% [fix_times, eye_xy] = find_events(text_file,'fixation');
% [sacc_end_times, sacc_end_xy, sacc_start_times, sacc_start_xy] = 
%    find_events(text_file,'saccade');
% message_time_and_number = find_events(text_file,'writeioport');
%
% INPUTS:
% - text_file should be the filename of the eyelink text file.
%    To create the eyelink text file, run the .edf file created during a
%    3D Search (3DS) experiment through the program Visual EDF2ASC.  For
%    event_type 'eyesample', use the EDF2ASC option 'samples only' or
%    'samples and events.'  Everything else requires 'events only' or
%    'samples and events.'
% - event_type is a string specifying the type of event you're looking for.
%
% EVENT TYPES AND CORRESPONDING OUTPUTS:
% - 'blink': output is the times at which a blink occurred.
% - 'button': output is the times at which a button was pressed.
% - 'eyesample': output 1 is the x and y position of the eye at each sample
%    output 2 is the size of the pupil at each sample
% - 'fixation': output 1 is the start and end time of each fixation
%    output 2 is the average x and y position of the eye during the fixation.
% - 'saccade': output 1 is the end time of the saccade
%    output 2 is the x and y position of the eye at the end of the saccade.
%    output 3 is the start time of the saccade, and output 4 is the x and y
%    position of the eye at the start of the saccade.
% - 'writeioport': output is an nx2 matrix where each row is [timestamp, 
%    number of message sent]. See GetNumbers.m for details.
%
% Created 11/12/10 by DJ.
% Updated 2/25/11 by DJ - added saccade start output
% Updated 5/31/11 by DJ - comments

% Set up
fid = fopen(text_file);
fseek(fid,0,'eof'); % find end of file
eof = ftell(fid);
fseek(fid,0,'bof'); % rewind to beginning

% Set the parameters for our search
% word: the word we check for to find relevant events
% format: the format we use in sscanf to turn a line of text into values of interest
% values: the values returned by sscanf (we specify the # of columns here)

% preallocate matrix size to speed up the loading process.
% For Jun.13th's experiment, ITI is accidentally set to 3.7-4.7s, the total
% max samples can hence be estimated as 1000(sampling rate)*125(trials)*4.7
% which arrives at 587500. we then set case 'eyesample' matrix size to
% zeros(587500,4).

switch event_type
    case 'blink'
        word = 'EBLINK'; 
        format = 'EBLINK %*s %d %d'; % Message format: EBLINK <eye (R/L)> <blinkstart> <blinkend>
        values = zeros(0,2);
    case 'button'
        word = 'BUTTON';
        format = 'BUTTON %d %*d %d'; % Message format: BUTTON <time> <button #> <state>
        values = zeros(0,2);
    case 'eyesample'
        word = '...';
        format = '%d %f %f %f'; % Message format: <time> <xpos> <ypos> <pupilsize> ...
        values = zeros(587500,4);    % do not skip the first time info
    case 'fixation'
        word = 'EFIX';
        format = 'EFIX %*s %d %d %*d %f %f %*f'; % Message format: EFIX <eye> <starttime> <endtime> <duration> <avgxpos> <avgypos> <avgpupilsize>
        values = zeros(0,4);
    case 'saccade'
        word = 'ESACC';
        format = 'ESACC %*s %d %d %*d %f %f %f %f'; % Message format: ESACC <eye> <starttime> <endtime> <duration> <startx> <starty> <endx> <endy> <amp> <peakvelocity>    
        values = zeros(0,6);
    case 'writeioport'
        word = 'write_ioport';
        format = 'MSG %d !CMD %*d write_ioport 0x378 %d'; % Message format: MSG <time> !CMD <delay???> write_ioport 0x378 <msgnumber>
        values = zeros(0,2);
        
    % record the messsages
    case 'synctime'
        word = 'SYNCTIME';
        format = 'MSG %d %*d SYNCTIME';
        values = zeros(0,1);
    case 'timeout'
        word = 'timeout';
        format = 'MSG %d %*d timeout';
        values = zeros(0,1);
    case 'image'
        word = '.bmp';
        format = 'MSG %d Trials %d/%d %s';
        values = zeros(0,2);      
    case 'sound'
        word = '.wav';
        format = 'MSG %d Trials %d/%d %s';
        values = zeros(0,2);
    otherwise
        error('FindEvents:InputType','event_type input %s not recognized!',event_type);
end


% let values update from the first row.
i=1;

% Get the messages we're looking for
% each row of 'values' will be the info for one line (e.g., timestamp, event)
while ftell(fid) < eof % if we haven't reached the end of the text file
    str = fgetl(fid); % read in next line of text file
    if findstr(str,word) % check for the code-word indicating a message was written
        if strcmp(word,'.bmp')==0 && strcmp(word,'.wav')==0
            stuff = sscanf(str,format)';
            if size(stuff,2)==size(values,2)
                values(i,:) = stuff;
                i = i+1;    % values = [values; stuff]; % add the info from this line as an additional row
            elseif strcmp(event_type,'eyesample')
                values(i,:) = [NaN,NaN,NaN,NaN];
                i = i+1;    % values = [values; NaN,NaN,NaN,NaN]; % add a blank sample so the time points still line up
            else
                warning('FindEvents:IncompleteEvent','Eyelink was unable to log the following event fully:\n %s',str); % sometimes saccades are not logged fully
            end
        elseif strcmp(word,'.bmp')==1
            stuff = sscanf(str,format)';
            % one of the message would return as 
            % MSG	758653 !V TRIAL_VAR cross greencircle.bmp, which gives
            % a stuff vector with only 1*1 size. Hence we need to filter
            % this out and set criterion for length for events.
            if length(stuff)>1 && stuff(4)==103 % the 4th element is 'g', which refers to 103.
                values = [values; [stuff(1), 2]]; % green = 2
            elseif length(stuff)>1 && stuff(4)==114 % the 4th element is 'r', which refers to 114.
                values = [values; [stuff(1), 1]]; % red = 1
            else
                warning('FindEvents:IncompleteEvent','Eyelink was unable to log the following event fully:\n %s',str); % sometimes saccades are not logged fully
            end
        elseif strcmp(word,'.wav')==1
            % auditory oddball doesn't include the above end message, hence
            % no need to worry.
            stuff = sscanf(str,format)';
            if stuff(4)==115 % the 4th element is 's', which refers to 115
                values = [values; [stuff(1), 2]]; % standard = 2
            elseif stuff(4)==118 % the 4th element is 'v', which refers to 118.
                values = [values; [stuff(1), 1]]; % oddball = 1
            else
                warning('FindEvents:IncompleteEvent','Eyelink was unable to log the following event fully:\n %s',str); % sometimes saccades are not logged fully
            end
        end            
    end
end

% Clean up
fclose(fid);



% Did we find any messages?
if size(values,1)==0 % if no
    warning('FindEvents:NoEventsFound',...
        'Couldn''t find ''%s'' anywhere in document ''%s''!',word,text_file);
    varargout = cell(1,nargout);
else
    % Crop the values into the desired output
    switch event_type
        case 'blink'
            varargout{1} = values(:,1); % output is times at which a blink started
            varargout{2} = values(:,2);
        case 'button'
            varargout{1} = values(values(:,2) == 1,1); % output is times at which any button was pressed
        case 'eyesample'
            varargout{1} = values(:,1);   % time of sample
            varargout{2} = values(:,2:3); % first output is the x and y position of the eye
            varargout{3} = values(:,4); % second output is the pupil size
        case 'fixation'
            varargout{1} = values(:,1:2); % first output is timestamp of start and end of fixation
            varargout{2} = values(:,3:4); % second output is avg. x and y position of eye during fixation
        case 'saccade'
            varargout{1} = values(:,2); % first output is timestamp of end of saccade
            varargout{2} = values(:,5:6); % second output is x and y position of eye at end of saccade
            varargout{3} = values(:,1); % third output is timestamp of start of saccade
            varargout{4} = values(:,3:4); % third output is timestamp of start of saccade
        case 'writeioport'
            varargout{1} = values; % output is [timestamp sent, message sent]
        
        case 'synctime'
            varargout{1} = values;
        case 'timeout'
            varargout{1} = values;
        case 'image'
            varargout{1} = values(:,1);
            varargout{2} = values(:,2);
        case 'sound'
            varargout{1} = values(:,1);
            varargout{2} = values(:,2);           
        otherwise 
            error('event_type input %s not recognized!',event_type);
    end
end