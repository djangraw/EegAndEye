function blinkRange = BlinkRange(x)

% Find the time ranges containing blinks in a 3DSarch or Squares data
% struct or EEGLAB file.
%
% blinkRange = BlinkRange(x)
% blinkRange = BlinkRange(EEG)
% 
% INPUTS:
% - x is a 3DSearch or Squares data struct.
% - EEG is an EEGLAB data struct for a 3DSearch or Squares data file.
%
% OUTPUTS:
% - blinkRange is an nx2 matrix, where n is the number of blinks in x.  The
% first column is the start time of each blink, and the second column is
% the end time of each blink.
%
% Created 8/8/11 by DJ.
% Updated 8/9/11 by DJ - eye blinks on the edges
% Updated 8/11/11 by DJ - allow EEG input, exclude blinks that are too long
% Updated 10/28/11 by DJ - make work for squares files
% Updated 11/28/11 by DJ - added text labeled events support

% define constants
MAX_BLINK_DURATION = 3000; % max allowable blink duration, in ms

% detect input type
if isfield(x,'setname')
    filetype = 'EEG';
    EEG = x;
else
    filetype = 'behavior';
end

% Extract Event Info
switch filetype
    case 'behavior'
        % Set up
        if isfield(x,'eyelink') % 3DSearch
            start_times = x.eyelink.saccade_start_times;
            end_times = x.eyelink.saccade_times;
            blink_times = x.eyelink.blink_times;
        elseif isfield(x,'target_color') % Squares
            start_times = x.saccade.start_time;
            end_times = x.saccade.end_time;
            blink_times = x.blink.start_time;
        end
        
    case 'EEG'
        % Only deal with continuous datasets
        if EEG.trials<1
            error('EEG dataset must be continuous!');
        end
        % extract event info
        eventTypes = [str2double({EEG.event(:).type})];
        eventLatencies = [EEG.event(:).latency] * 1000/EEG.srate;
        
        % Detect file type
        if strcmp(EEG.setname(1:3),'3DS')
            GetNumbers;
            start_times = eventLatencies(eventTypes==Numbers.SACCADE_START);
            end_times = eventLatencies(eventTypes==Numbers.SACCADE_END);
            blink_times = eventLatencies(eventTypes==Numbers.BLINK);           
        elseif strcmp(EEG.setname(1:2),'sq')
            Constants = GetSquaresConstants;
            if mean(isnan(eventTypes))<0.5 % Mostly numeric codes
                start_times = eventLatencies(ismember(eventTypes,Constants.SACCADESTART_BASE:Constants.SACCADEEND_BASE));
                end_times = eventLatencies(ismember(eventTypes,Constants.SACCADEEND_BASE+(0:(Constants.SACCADEEND_BASE-Constants.SACCADESTART_BASE))));
                blink_times = eventLatencies(eventTypes==Constants.BLINKSTART);   
            else % String codes
                eventTypes = {EEG.event(:).type};
                start_times = eventLatencies(strmatch('SS',eventTypes));
                end_times = eventLatencies(strmatch('SE',eventTypes));
                blink_times = eventLatencies(strmatch('BlinkStart',eventTypes));  
            end
        else
            error('File type not detected!')
        end
end
        

% get rid of eye blinks on the edges
blink_times = blink_times(blink_times>start_times(1) & blink_times<end_times(end));

% Get ranges
nBlinks = length(blink_times);
blinkRange = nan(nBlinks,2);
for i=1:nBlinks    
    blinkRange(i,1) = start_times(find(start_times<blink_times(i),1,'last')); % last saccade start event before blink
    blinkRange(i,2) = end_times(find(end_times>blink_times(i),1,'first')); % first saccade end event after blink       
end

% get rid of blinks where the saccade around it was not properly logged
% (e.g. on the edge of concatenated sessions)
blinkRange = blinkRange(blinkRange(:,2)-blinkRange(:,1)<MAX_BLINK_DURATION,:); % exclude blinks that are too long to be real