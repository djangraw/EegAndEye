% 3D Search Data Extraction Tool
% Run this script to load the three data files that are output from a 
% 3dsearch experiment and compare the data they contain.
% Checks for inconsistencies in timing and outputs a data struct x with 
% everything(?) in it.
%
% Created 6/14/10 by DJ.
% Updated 6/17/10 by DJ - added lots of stuff.


%% GET FILENAMES
eye_filename = uigetfile('*.txt', 'Select Eye Tracker Text File');
eeg_filename = uigetfile('*.bdf', 'Select EEG File');

%% READ UNITY LOG PARAMETERS
disp('Getting Unity Session Parameters...')

% titles from text file
params = {'EDF_filename:', 'subject:', 'session:', 'level:', ...
    'trialTime:', 'targetFolder:', 'distractorFolder:', 'locations:',...
    'objectSize:', 'distanceToPopup:', 'objectPrevalence:',... 
    'targetPrevalence:','minPopupDelay:', 'maxPopupDelay:',...
    'objectMoveTime:', 'popupStayTime:', 'GuiIsOn:', 'pixAround:',...
    'eyelink.offset_x:', 'eyelink.offset_y:', 'eyelink.gain_x:',...
    'eyelink.gain_y:'};
% corresponding field name in struct we will create
fields = {'EDF_filename', 'subject', 'session', 'level', ...
    'trial_time', 'target_category', 'distractor_category', 'locations',...
    'object_size', 'distance_to_popup', 'object_prevalence',...
    'target_prevalence', 'min_popup_delay', 'max_popup_delay',...
    'object_move_time', 'popup_stay_time', 'gui_is_on','pix_around',...
    'eye_offset_x', 'eye_offset_y', 'eye_gain_x', 'eye_gain_y'};

% Take each tag from "params", find it in the Unity text file, and move the
% value that follows it into the corresponding field in the struct x.
x = struct(); % create empty struct
for i=1:numel(params)
    results = find_word(eye_filename,params{i},0,'%s'); % get the value following the tag
    result = results{1}; % grab the first one
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

%% READ EYE TRACKER FILE
% This is used for events like objects entering/exiting the scene, 
% saccades, and for checking the eyelink and eeg times against each other.
disp('Getting Eye Tracker Events...');

% Get Port Events
[port_events starttime] = find_writeioport(eye_filename);

% Get Saccade Times
saccade_times = find_saccadeevents(eye_filename);

% Get Saccade Events
saccade_events = matchSaccadesToObjects(saccade_times, port_events);

% Put results in struct
x.eyelink.starttime = starttime;
x.eyelink.port_events = port_events;
x.eyelink.saccade_events = saccade_events;

% Clean up
clear port_events starttime saccade_times saccade_events
disp('...Success!')

%% GET TRIAL INFO
disp('Getting trial info from eyelink data file...')
[trialtimes trialtypes] = get_trialtimes(x.eyelink.port_events);
nTrials = size(trialtimes,1);

trialobjects = get_trialobjects(eye_filename);

if numel(trialobjects) ~= nTrials 
    error('Trial output for times and objects do not match!');
end

for i=1:nTrials
    x.trials(i).loadtime = trialtimes(i,1);
    x.trials(i).starttime = trialtimes(i,2);
    x.trials(i).endtime = trialtimes(i,3);
    x.trials(i).type = trialtypes(i); % surpriseLevel for the trial
    x.trials(i).objects = trialobjects{i};
end

% Clean up
clear trialtimes nTrials trialobjects trialtypes i
disp('...Success!')

%% GET EEG EVENTS
disp('Getting EEG Parallel Port Events...')

x.eeg.port_events = get_biosemi_events(eeg_filename);
x.eeg.eventsamplerate = get_biosemi_eventsamplerate(eeg_filename);
disp('...Success!')



%% COMPARE EYE AND EEG EVENTS
disp('Checking Parallel Port Events (Eye Tracker vs. EEG)...');
GetNumbers; % Load the 'Numbers' struct, which has the latest Unity constants

% Extract info from the struct we made (for easier access)
eye = x.eyelink.port_events; % should start from the 'START_RECORDING' broadcast
eeg = x.eeg.port_events(find(x.eeg.port_events(:,2)==Numbers.START_RECORDING,1):end,:); % should start from the 'START_RECORDING' broadcast

% Take care of sampling rates, offsets
eeg(:,1) = eeg(:,1)*1000/x.eeg.eventsamplerate; % to get time in ms instead of samples (EEG sampling rate = 2048)
eeg(:,1) = eeg(:,1) - eeg(1,1); % make time of first event t=0
eye(:,1) = eye(:,1) - eye(1,1); % make time of first event t=0

% Make sure events are the same
if size(eeg,1) ~= size(eye,1) || sum(eeg(:,2)~=eye(:,2)) > 0 % number of events that aren't in exactly the right spot
    disp('Events don''t match up!');
else
    disp('All event types match. Checking timing...')
    delay = eeg(:,1)-eye(:,1); % subtract timestamps
    figure(1); clf;
    subplot(2,1,1)
    hist(delay);
    xlabel('Delay (Time Received - Time Sent) in ms')
    ylabel('# events')
    title('Event Timing Consistency Check');
    subplot(2,1,2);
    plot(eye(:,1)/1000,delay, '.');
    xlabel('Event time (s)');
    ylabel('Delay in ms');
    disp(sprintf('average absolute value of delay is %g', mean(abs(delay))));
end

clear wrong_type delay eye eeg Numbers


%% VIEW LIFETIMES OF OBJECTS

for i=1:numel(x.trials)
    PlotObjectLifetimes(x.trials(i),x.eyelink.saccade_events,x.eyelink.port_events)
end



