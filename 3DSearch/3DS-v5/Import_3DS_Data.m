function x = Import_3DS_Data(subject,session)

% x = Import_3DS_Data(subject,session)
%
% 3D Search Data Extraction Tool
% Run this function to load the three data files that are output from a 
% 3dsearch experiment and compile the data they contain into a struct.
% Checks for inconsistencies in timing and outputs a data struct x with 
% everything in it.  Saves the result in the current directory as
% '3DS-[subject]-[session].mat'.
%
% Created 6/14/10 by DJ (called CheckData)
% Updated 6/17/10 by DJ - added lots of stuff (called CheckData)
% Updated 7/22/10 by DJ - adjusted to fit new output format, made it a fn.
% Updated 7/27/10 by DJ - changed field names and many helper functions
% Updated 7/29/10 by DJ - changed events field back to eyelink, moved
%   trial_type, added eyelink to eeg time conversion.

tic

%% GET FILENAMES
disp(sprintf('---Importing Data For Subject %d, Session %d---',subject,session));
eye_filename = sprintf('3DS-%d-%d.asc',subject,session);
event_filename = sprintf('post-3DS-%d-%d.txt',subject,session);
eeg_filename = sprintf('3DS-%d-%d.bdf',subject,session);
    
%% READ EYELINK LOG PARAMETERS
disp('Getting Real-Time Session Parameters...')

% titles from text file
params = {'EDF_filename:', 'subject:', 'session:', 'level:', ...
    'SurpriseLevel:','trialTime:', 'targetFolder:', 'distractorFolder:',...
    'locations:','objectSize:', 'distanceToPopup:', 'objectPrevalence:',... 
    'targetPrevalence:','minPopupDelay:', 'maxPopupDelay:',...
    'objectMoveTime:', 'popupStayTime:', 'GuiIsOn:'};
% corresponding field name in struct we will create
fields = {'EDF_filename', 'subject', 'session', 'level', ...
    'trial_type','trial_time', 'target_category', 'distractor_category',...
    'locations','object_size', 'distance_to_popup', 'object_prevalence',...
    'target_prevalence', 'min_popup_delay', 'max_popup_delay',...
    'object_move_time', 'popup_stay_time', 'gui_is_on'};

% Take each tag from "params", find it in the Unity text file, and move the
% value that follows it into the corresponding field in the struct x.
x = struct(); % create empty struct
for i=1:numel(params)
    result = find_first_word(eye_filename,params{i},'%s'); % get the value following the tag
    numresult = str2num(result); % is it a number or a string?
    if isempty(numresult) % if it's a string
        x = setfield(x,fields{i},result);
    else % if it's a number
        x = setfield(x,fields{i},numresult);
    end
        
end

% Clean up
clear fields params results result numresult i
disp('...Success!')

%% READ EVENT LOG PARAMETERS
disp('Getting Post-Hoc Session Parameters...')

% titles from text file
params = {'eyelink.offset_x:', 'eyelink.offset_y:', ...
    'eyelink.gain_x:', 'eyelink.gain_y:'};
% corresponding field name in struct we will create
fields = {'eye_offset_x', 'eye_offset_y', ...
    'eye_gain_x', 'eye_gain_y'};

% Take each tag from "params", find it in the Unity text file, and move the
% value that follows it into the corresponding field in the struct x.
x.calibration = struct(); % create empty struct
for i=1:numel(params)
    result = find_first_word(event_filename,params{i},'%s'); % get the value following the tag
    numresult = str2num(result); % is it a number or a string?
    if isempty(numresult) % if it's a string
        x.calibration = setfield(x.calibration,fields{i},result);
    else % if it's a number
        x.calibration = setfield(x.calibration,fields{i},numresult);
    end
        
end

% Clean up
clear fields params results result numresult i
disp('...Success!')

%% READ EVENTS FILE
% This is used for events like objects entering/exiting the scene, 
% saccades to objects, etc.
disp('Getting events file info...')
object_events = get_objectevents(event_filename);

% Put results in struct
x.eyelink.object_events = object_events;

% Clean up
clear object_events
disp('Success!')

%% READ EYE TRACKER FILE
% This is used for events like objects entering/exiting the scene, 
% saccades, and for checking the eyelink and eeg times against each other.
disp('Getting Eye Tracker Events...');

% Get Port Events
sync_events = find_writeioport(eye_filename);

button_times = find_buttonevents(eye_filename);

% Put results in struct
x.eyelink.sync_events = sync_events;
x.eyelink.button_times = button_times;

% Clean up
clear sync_events button_times
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
x.eyelink.load_time = trialtimes(1);
x.eyelink.start_time = trialtimes(2);
x.eyelink.end_time = trialtimes(3);
x.eyelink.object_lifetimes = lifetimes{1};
x.objects = objects{1};

% Clean up
clear trialtimes nTrials trialobjects lifetimes
disp('...Success!')

%% GET EEG EVENTS
disp('Getting EEG Parallel Port Events...')

% Get event sample rate (typically 2048 Hz for Biosemi)
x.eeg.eventsamplerate = get_biosemi_eventsamplerate(eeg_filename);
% Get events
GetNumbers;
sync_events = get_biosemi_events(eeg_filename); % Get raw events
x.eeg.sync_events = sync_events(find(sync_events(:,2)==Numbers.START_RECORDING,1):end,:); % should start from the 'START_RECORDING' broadcast

disp('...Success!');
clear sync_events Numbers;

%% CONVERT EYELINK EVENTS TO EEG-BASED TIMES

disp('Converting eyelink events to eeg times..');

x = EyelinkToEegTimes(x);

disp('...Success!');

%% COMPARE EYE AND EEG SYNC EVENTS
disp('Checking Parallel Port Events (Eye Tracker vs. EEG)...');

figure(session+20); clf;
CheckSyncEvents(x);

disp('Success!');

%% VIEW LIFETIMES OF OBJECTS

figure(session); clf;
PlotObjectLifetimes(x) % as a check that we've imported our data properly

clear i;

%% SAVE STRUCT
disp('Saving result in current directory...')
save(sprintf('3DS-%d-%d.mat',subject,session),'x');
disp('...Success!')
toc