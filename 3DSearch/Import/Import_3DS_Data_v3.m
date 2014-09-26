function x = Import_3DS_Data_v3(subject,ascSession,eegSession,eegFiletype,pixelThresholds,timeLimits)

% x = Import_3DS_Data_v3(subject,ascSession,eegSession,eegFiletype,pixelThresholds,timeLimits)
%
% 3D Search Data Extraction Tool
% Run this function to load the three data files that are output from a 
% 3dsearch experiment and compile the data they contain into a struct.
% Checks for inconsistencies in timing and outputs a data struct x with 
% everything in it.  Saves the result in the current directory as
% '3DS-[subject]-[ascSession].mat'.
%
% INPUTS:
% -
%
% OUTPUTS:
% -x is the struct with all the info in it.
%
% Created 6/14/10 by DJ (called CheckData)
% Updated 6/17/10 by DJ - added lots of stuff (called CheckData)
% Updated 7/22/10 by DJ - adjusted to fit new output format, made it a fn.
% Updated 7/27/10 by DJ - changed field names and many helper functions
% Updated 7/29/10 by DJ - changed events field back to eyelink, moved
%   trial_type, added eyelink to eeg time conversion.
% Updated 8/18/10 by DJ - changed eye filename to .asc, added
%   object_limits, saccade_times, and saccade_position to reflect new
%   output of Unity replays (3DSearch_5pt3). Commented.
% Updated 8/25/10 by DJ - changed str2num to str2double, added dynamic
%   struct fields, and used fprintf.
% Updated 8/31/10 by DJ - changed param_field structure, and added fields.
% Updated 9/13/10 by DJ - added fixation info
% Updated 9/21/10 by DJ - added active_session, move_speed, and spin_speed,
%   switched trial_type to surprise_level.
% Updated 10/18/10 by DJ - added blink times.
% Updated 11/16/10 by DJ - switched all find_xxx to find_events.
% Updated 11/23/10 by DJ - added filetype switch to accommodate sensorium
% Updated 12/15/10 by DJ - added ascSession/eegSession
% Updated 2/23/11 by DJ - added inputs pixelThresholds and timeLimts. TO DO: update comments!
% Updated 2/25/11 by DJ - added saccade_start fields
% Updated 7/28/11 by DJ - added new Sensorium cap option
% (eegFileType='Sensorium-2011'), added brake button presses
% Updated 8/15/11 by DJ - made button events date-specific
% Updated 10/26/11 by DJ - cropped sync events to Numbers.SYNC,0 only!

%% SET UP
tic
if nargin<3 || strcmp(eegFiletype,'');
    eegFiletype = 'Biosemi';
end
if nargin<4 || isempty(pixelThresholds)
    pixelThresholds = 1000; % include all saccades
end
if nargin<5 || isempty(timeLimits)
    timeLimits = [0 Inf]; % include all saccades
end

%% GET FILENAMES
fprintf('---Importing Data For Subject %d, Session %d---\n',subject,ascSession);
eye_filename = sprintf('3DS-%d-%d.asc',subject,ascSession); % eyelink .edf file converted to text using Visual EDF2ASC.
event_filename = sprintf('post-3DS-%d-%d.txt',subject,ascSession); % text file written during a Unity replay
switch eegFiletype
    case 'Biosemi'
        eeg_filename = sprintf('3DS-%d-%d.bdf',subject,eegSession); % .bdf EEG file.
    case {'Sensorium-2008' 'Sensorium-2011' 'Sensorium-2011-50Hz' 'Sensorium-2011-1000Hz'}
        eeg_filename = sprintf('3DS-%d.%.3d.b',subject,eegSession); % .dat EEG file
end
    
%% READ EYELINK LOG PARAMETERS
disp('Getting Real-Time Session Parameters...')

% Translate from headings in Eyelink messages to fields in matlab struct.
% Each row is, for one field:
% Titles from text file, corresponding field name in struct we'll create
param_field = {'EDF_filename:', 'EDF_filename';...
    'subject:', 'subject';...
    'session:', 'session';...
    'Date:', 'date';...
    'replayFilename:', 'replayFilename';...    
    'level:', 'level';...
    'SurpriseLevel:', 'surprise_level';...
    'isActiveSession:', 'active_session';...
    'trialTime:', 'trial_time';...
    'targetFolder:', 'target_category';...
    'distractorFolder:', 'distractor_category';...
    'locations:', 'locations';...
    'objectSize:', 'object_size';...
    'distanceToPopup:', 'distance_to_popup';...
    'objectPrevalence:', 'object_prevalence';...
    'targetPrevalence:', 'target_prevalence';...
    'minPopupDelay:', 'min_popup_delay';...
    'maxPopupDelay:', 'max_popup_delay';...
    'objectMoveTime:', 'object_move_time';...
    'popupStayTime:', 'popup_stay_time';...
    'startPoint:', 'start_point';...
    'moveSpeed:', 'move_speed';...
    'spinSpeed:', 'spin_speed'};
 

% Look in the "SESSION PARAMETERS" section only
[~, endpoint] = find_first_word(eye_filename,'END'); % find the point where the session parameters end

% Take each tag from "params", find it in the Unity text file, and move the
% value that follows it into the corresponding field in the struct x.
x = struct(); % create empty struct
for i=1:size(param_field,1)
    result = find_first_word(eye_filename,param_field{i,1},'%s',endpoint); % get the value following the tag
    numresult = str2double(result); % is it a number or a string?
    if isnan(numresult) % if it's a string
        x.(param_field{i,2}) = result;
    else % if it's a number
        x.(param_field{i,2}) = numresult;
    end
        
end

% Clean up
clear param_field results result numresult endpoint i
disp('...Success!')

%% READ EVENT LOG PARAMETERS
disp('Getting Post-Hoc Session Parameters...')

param_field = {'eyelink.offset_x:', 'eye_offset_x';...
    'eyelink.offset_y:', 'eye_offset_y';...
    'eyelink.gain_x:', 'eye_gain_x';...
    'eyelink.gain_y:', 'eye_gain_y'};

% Look in the "SESSION PARAMETERS" section only
[~, endpoint] = find_first_word(event_filename,'LOAD'); % find the point where the session parameters end


% Take each tag from "params", find it in the Unity text file, and move the
% value that follows it into the corresponding field in the struct x.
x.calibration = struct(); % create empty struct
for i=1:size(param_field,1)
    result = find_first_word(event_filename,param_field{i,1},'%s',endpoint); % get the value following the tag
    numresult = str2double(result); % is it a number or a string?
    if isnan(numresult) % if it's a string
        x.calibration.(param_field{i,2}) = result;
    else % if it's a number
        x.calibration.(param_field{i,2}) = numresult;
    end
        
end

% Clean up
clear param_field results result numresult endpoint i
disp('...Success!')

%% READ EVENTS FILE
% This is used for events like objects entering/exiting the scene, 
% saccades to objects, etc.
disp('Getting events file info...')
object_limits = get_objectlimits(event_filename);
object_events = get_objecttimes(object_limits);
% object_events = get_objectevents(event_filename);

% Put results in struct
x.eyelink.object_limits = object_limits; % [time object left top width height] for each frame where an object was visible.
x.eyelink.object_events = object_events; % timestamp and event number for object entrances and exits

% Clean up
clear object_limits object_events pixAround
disp('...Success!')

%% READ EYE TRACKER FILE
% This is used for events like objects entering/exiting the scene, 
% saccades, and for checking the eyelink and eeg times against each other.
disp('Getting Eye Tracker Events...')

GetNumbers(x.date);
% Get Port Events
sync_events = find_events(eye_filename, 'writeioport');
% Get button presses
[all_button_times all_button_numbers] = find_events(eye_filename, 'button');
button_times = all_button_times(ismember(all_button_numbers,Numbers.GAMEPAD_TARGET_BUTTON));
brake_button_times = all_button_times(ismember(all_button_numbers,Numbers.GAMEPAD_BRAKE_BUTTON));
% Get saccade and fixation info
[saccade_times endpos saccade_start_times startpos] = find_events(eye_filename, 'saccade');
[fixation_times avgpos] = find_events(eye_filename, 'fixation');
% Convert eye position from eyelink coordinates to Unity coordinates
% endpos(:,1) = (endpos(:,1) - x.calibration.eye_offset_x) * x.calibration.eye_gain_x;
% endpos(:,2) = (endpos(:,2) - x.calibration.eye_offset_y) * x.calibration.eye_gain_y;
% startpos(:,1) = (startpos(:,1) - x.calibration.eye_offset_x) * x.calibration.eye_gain_x;
% startpos(:,2) = (startpos(:,2) - x.calibration.eye_offset_y) * x.calibration.eye_gain_y;
% avgpos(:,1) = (avgpos(:,1) - x.calibration.eye_offset_x) * x.calibration.eye_gain_x;
% avgpos(:,2) = (avgpos(:,2) - x.calibration.eye_offset_y) * x.calibration.eye_gain_y;
endpos = ApplyEyeCalibration(endpos,x);
startpos = ApplyEyeCalibration(startpos,x);
avgpos = ApplyEyeCalibration(avgpos,x);
% Get blink times
blink_times = find_events(eye_filename, 'blink');
% crop sync events
sync_events = sync_events(find(sync_events(:,2)==Numbers.START_RECORDING,1):end,:); % should start from the 'START_RECORDING' broadcast
sync_events = sync_events(ismember(sync_events(:,2),[Numbers.START_RECORDING, Numbers.SYNC,0]),:);

% Put results in struct
x.eyelink.sync_events = sync_events; % timestamp and event number of sync events
x.eyelink.button_times = button_times; % timestamp of target button presses
x.eyelink.brake_button_times = brake_button_times; % timestamp of brake button presses
x.eyelink.saccade_times = saccade_times; % timestamp of saccade endpoints
x.eyelink.saccade_positions = endpos; % x and y position (in Unity coordinates) of saccade endpoints
x.eyelink.saccade_start_times = saccade_start_times; % timestamp of saccade endpoints
x.eyelink.saccade_start_positions = startpos; % x and y position (in Unity coordinates) of saccade endpoints
x.eyelink.fixation_times = fixation_times; % timestamp of fixation start and end
x.eyelink.fixation_positions = avgpos; % avg. x and y position (in Unity coordinates) of fixation
x.eyelink.blink_times = blink_times; % timestamp of eye blinks

% Clean up
clear sync_events *button_* saccade_times endpos fixation_times avgpos blink_times Numbers
disp('...Success!')

%% GET OBJECT EVENTS
disp('Getting Saccade Events...')
saccade_events = classify_saccades(x.eyelink.saccade_times,x.eyelink.saccade_positions,x.eyelink.object_limits,pixelThresholds,timeLimits);

x.eyelink.saccade_events = saccade_events;

clear saccade_events
disp('...Success!')

%% GET TRIAL INFO
disp('Getting trial info from eyelink data file...')
trialtimes = get_trialtimes(eye_filename);
nTrials = size(trialtimes,1);
if nTrials>1
    error('More than one trial in given session!');
end

[objects lifetimes] = get_trialobjects(eye_filename);

if numel(objects) ~= nTrials 
    error('Trial output for times and objects do not match!');
end

% Load info into struct
x.eyelink.record_time = find_first_word(eye_filename,'START','%d'); % timestamp of first events recorded by Eyelink
x.eyelink.load_time = trialtimes(1); % timestamp of when the trial and its objects started being loaded
x.eyelink.start_time = trialtimes(2); % timestamp of when loading was complete and the subject could start moving
x.eyelink.end_time = trialtimes(3); % timestamp of when the trial ended and the display froze
x.eyelink.object_lifetimes = lifetimes{1}; % creation and destruction times of each object
x.objects = objects{1};

% Clean up
clear trialtimes nTrials trialobjects lifetimes
disp('...Success!')

%% GET EEG EVENTS

switch eegFiletype
    case 'Biosemi'
        % Get event sample rate (typically 2048 Hz for Biosemi)
        x.eeg.eventsamplerate = get_biosemi_eventsamplerate(eeg_filename);
        % Get events
        GetNumbers;
        sync_events = get_biosemi_events(eeg_filename); % Get raw events
        x.eeg.sync_events = sync_events(find(sync_events(:,2)==Numbers.START_RECORDING,1):end,:); % should start from the 'START_RECORDING' broadcast
    case {'Sensorium-2008' 'Sensorium-2011'}
        x.eeg.eventsamplerate = 1000;
        % Get events
        GetNumbers;
        sync_events = get_sensorium_events(eeg_filename); % Get raw events
        if subject==13 && ascSession==4 % TEMP FIX FOR AN INCOMPLETE DATASET!
            sync_events = [-869 200; sync_events];
        elseif subject==13 && ascSession==9 % TEMP FIX FOR AN INCOMPLETE DATASET!
            x.eyelink.sync_events = x.eyelink.sync_events(1:size(sync_events,1),:);
        end
        sync_events = sync_events(find(sync_events(:,2)==Numbers.START_RECORDING,1):end,:); % should start from the 'START_RECORDING' broadcast
        x.eeg.sync_events = sync_events(ismember(sync_events(:,2),[Numbers.START_RECORDING, Numbers.SYNC,0]),:);
end
        
disp('...Success!');
clear sync_events Numbers

%% CONVERT EYELINK EVENTS TO EEG-BASED TIMES

disp('Converting eyelink events to eeg times..');

x = EyelinkToEegTimes(x);
disp('...Success!');

%% COMPARE EYE AND EEG SYNC EVENTS
disp('Checking Parallel Port Events (Eye Tracker vs. EEG)...');

figure(ascSession+20); clf;
CheckSyncEvents(x);

disp('Success!');

%% VIEW LIFETIMES OF OBJECTS

figure(ascSession+40); clf;
PlotObjectLifetimes(x) % as a check that we've imported our data properly

clear i;

%% SAVE STRUCT
disp('Saving result in current directory...')
save(sprintf('3DS-%d-%d.mat',subject,ascSession),'x');
disp('...Success!')
toc % Display elapsed time