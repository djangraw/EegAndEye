function varargout = find_events(text_file,event_type)

% [varargout] = find_events(text_file,event_type)
%
% EXAMPLE USAGE:
% [blink_start_times, blink_end_times] = find_events(text_file,'blink');
% [button_press_times, button_number] = find_events(text_file,'button');
% [eye_xy, pupil_size] = find_events(text_file,'eyesample');
% [fix_times, eye_xy] = find_events(text_file,'fixation');
% [sacc_end_times, sacc_end_xy, sacc_start_times, sacc_start_xy] = 
%    find_events(text_file,'saccade');
% message_time_and_number = find_events(text_file,'writeioport');
% [message_time trial_type] = find_events(text_file,'trialtype');
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
% - 'blink': output is the times at which a blink began and ended.
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
% - 'trialtype': output 1 is the time when the trialtype message was sent
%    (after a trial ends in the Squares experiment)
%    output 2 is the trial type, a six-digit number indicating the trial 
%    type.  The first digit is 1 if the fixation cross was on the left and 
%    2 if it was on the right.  Digits 2-6 indicate the color code of each 
%    square from right to left.  A code of 2 or 4 indicates the middle bar 
%    of that square was blue. See ExperimentBuilder code for details.
% - GENERIC display option: if it's not any of the above, the program will
%    search for 'MSG <x> <y> <event_type>', and will return x-y (the time
%    when eyelink reports that it displayed the event in question).
%
% Created 11/12/10 by DJ.
% Updated 2/25/11 by DJ - added saccade start output
% Updated 5/31/11 by DJ - comments
% Updated 7/28/11 by DJ - added 2nd output for button option
% Updated 10/17/11 by DJ - added trialtype option for Squares experiment
% Updated 3/4/13 by DJ - added generic display option ('otherwise' in code)
% Updated 12/5/13 by DJ - added leader option

% Set up
fid = fopen(text_file);
fseek(fid,0,'eof'); % find end of file
eof = ftell(fid);
fseek(fid,0,'bof'); % rewind to beginning

% Set the parameters for our search
% word: the word we check for to find relevant events
% format: the format we use in sscanf to turn a line of text into values of interest
% values: the values returned by sscanf (we specify the # of columns here)
switch event_type
    case 'blink'
        word = 'EBLINK'; 
        format = 'EBLINK %*s %d %d'; % Message format: EBLINK <eye (R/L)> <blinkstart> <blinkend>
        values = zeros(0,2);
    case 'button'
        word = 'BUTTON';
        format = 'BUTTON %d %d %d'; % Message format: BUTTON <time> <button #> <state>
        values = zeros(0,3);
    case 'eyesample'
        word = '...';
        format = '%*d %f %f %f'; % Message format: <time> <xpos> <ypos> <pupilsize> ...
        values = zeros(0,3);
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
    case 'trialtype'
        word = 'Trial';
        format = 'MSG %d Trial %*d Type %d';
        values = zeros(0,2);       
    case 'leader'
        word = 'Leader';
        format = 'MSG %d Leader %s';
        values = zeros(0,2);
    otherwise
        warning('FindEvents:InputType','event_type input %s not recognized!',event_type);
        word = event_type;
        format = ['MSG %d %d ' event_type];
        values = zeros(0,2);
end



% Get the messages we're looking for
% each row of 'values' will be the info for one line (e.g., timestamp, event)
while ftell(fid) < eof % if we haven't reached the end of the text file
    str = fgetl(fid); % read in next line of text file
    if findstr(str,word) % check for the code-word indicating a message was written
        stuff = sscanf(str,format)';
        if size(stuff,2)==size(values,2)
            values = [values; stuff]; % add the info from this line as an additional row
        elseif strcmp(event_type,'eyesample')
            values = [values; NaN,NaN,NaN]; % add a blank sample so the time points still line up
        elseif strcmp(event_type,'leader')
            if strcmpi(char(stuff(2:end)),'Slow')
                values = [values; stuff(1), 1]; % second item will be 'is this slow (1) or fast (2)?'
            elseif strcmpi(char(stuff(2:end)),'Fast')
                values = [values; stuff(1), 2]; 
            end
        else
            warning('FindEvents:IncompleteEvent','Eyelink was unable to log the following event fully:\n %s',str); % sometimes saccades are not logged fully
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
            varargout{1} = values(:,1); % first output is times at which a blink started
            varargout{2} = values(:,2); % second output is times at which a blink ended
        case 'button'
            varargout{1} = values(values(:,3) == 1,1); % first output is times at which any button was pressed
            varargout{2} = values(values(:,3) == 1,2); % second output is number of the button was pressed
        case 'eyesample'
            varargout{1} = values(:,1:2); % first output is the x and y position of the eye
            varargout{2} = values(:,3); % second output is the pupil size
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
        case 'trialtype'
            varargout{1} = values(:,1); % first output is Trial message time (sent after trial ends)
            varargout{2} = values(:,2); % second output is Trial type number        
        case 'leader'            
            varargout{1} = values(values(:,2)==1,1); % slow times
            varargout{2} = values(values(:,2)==2,1); % fast times
        otherwise 
%             error('event_type input %s not recognized!',event_type);
            varargout{1} = values(:,1)-values(:,2); % first output is display time
    end
end