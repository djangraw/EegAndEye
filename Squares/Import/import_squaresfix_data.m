function x = import_squaresfix_data(subject, ascSession, eegSession, eegFiletype, pixel_threshold, drift_correct,prefix)

% Import data from a squares experiment
%
% x = import_squaresfix_data(subject, ascSession, eegSession, pixel_threshold, drift_correct,prefix)
%
% INPUTS:
% - subject is a scalar indicating the subject number in the filenames
% - ascSession is a scalar indicating the session number of the eyelink
% (.asc) datafile.
% - eegSession is a scalar indicating the session number of the eeg (.bdf
% or .b.dat or .bh.dat) datafile.
% - pixel_threshold is a scalar indicating the distance in pixels between
% the center of a square and the mean fixation point that will still be 
% considered a fixation on that square.  (currently squares are 150 pixels
% apart).
% - drift_correct is a binary value indicating whether you want to offset
% the eye position on each trial by assuming the closest saccade during the
% fixation period was at the center of the fixation cross
% (see perform_drift_correction.m for details).
%
% OUTPUTS:
% - x is a struct with lots of info about each trial, saccade, etc. in it.
% See code for details.
%
% Created 3/4/13 by DJ based on import_squares_data.m.
% Updated 3/19/13 by DJ - added trialnum for saccades & fixations
% Updated 5/17/13 by DJ - added prefix input, sf3 support

% Set up
if ~exist('pixel_threshold','var')
    pixel_threshold = 50;
end
if ~exist('prefix','var')
    prefix = 'sf';
end

%% GET FILENAMES
fprintf('---Importing Data For Subject %d, Session %d---\n',subject,ascSession);
eye_filename = sprintf('%s_%d_%d.asc',prefix,subject,ascSession); % eyelink .edf file converted to text using Visual EDF2ASC.

switch eegFiletype
    case 'Biosemi'
        eeg_filename = sprintf('%s-%d-%d.bdf',prefix,subject,eegSession); % .bdf EEG file.
    case {'Sensorium-2008' 'Sensorium-2011' 'Sensorium-2011-1000Hz' 'Sensorium-2011-50Hz'}
        eeg_filename = sprintf('%s-%d.%.3d.b',prefix,subject,eegSession); % .dat EEG file
    case 'Sensorium-2013'
        eeg_filename = sprintf('%s-%d.%.3d.bh',prefix,subject,eegSession); % .dat EEG file
end

% Initialize struct with this info
x.subject = subject;
x.session = ascSession;
x.eyeFilename = eye_filename;
x.eegFilename = eeg_filename;
x.experiment = find_first_word(eye_filename,'BY','%s');

%% READ EYE TRACKER FILE
% This is used for events like objects entering/exiting the scene, 
% saccades, and for checking the eyelink and eeg times against each other.
disp('Getting Eye Tracker Events...')

% Get recording start time
recording_start_time = find_first_word(eye_filename,'START','%d');

% Get Port Events
Constants = GetSquaresConstants;
sync_events = find_events(eye_filename, 'writeioport');
block_start_time = sync_events(sync_events(:,2)==Constants.BLOCK_START,1);
trial_fix_time = sync_events(sync_events(:,2)==Constants.FIXCROSS_ON,1);
trial_start_time = sync_events(sync_events(:,2)==Constants.TRIAL_START,1);
trial_end_time = sync_events(sync_events(:,2)==Constants.TRIAL_END,1);
trial_circle_time = sync_events(sync_events(:,2)==Constants.ENDCIRCLE_OFF,1);
if isempty(trial_circle_time)
    trial_circle_time = trial_end_time;
end
% get trial types
[~,trial_type] = find_events(eye_filename, 'trialtype');

% Get button presses
[button_time button_number] = find_events(eye_filename, 'button');
% Get saccade and fixation info
[saccade_end_time endpos saccade_start_time startpos] = find_events(eye_filename, 'saccade');
[fixation_time avgpos] = find_events(eye_filename, 'fixation');
% Get blink times
[blink_start_time, blink_end_time] = find_events(eye_filename, 'blink');

% Get square display times
display_square_time = find_events(eye_filename, 'DISPLAY_SQUARE');

% Put results in struct
x.sync.eyelink = sync_events(:,1);
x.sync.events = sync_events(:,2); % timestamp and event number of sync events
x.recording_start_time = recording_start_time;
x.block_start_time = block_start_time;
x.trial.fix_time = trial_fix_time;
x.trial.start_time = trial_start_time;
x.trial.end_time = trial_end_time;
x.trial.circle_time = trial_circle_time;
x.trial.type = trial_type;
x.button.time = button_time; % timestamp of button presses
x.button.number = button_number; % number of button pressed
x.saccade.start_time = saccade_start_time; % timestamp of saccade startpoint
x.saccade.start_position = startpos; % x and y position of saccade startpoint
x.saccade.end_time = saccade_end_time; % timestamp of saccade endpoint
x.saccade.end_position = endpos; % x and y position of saccade endpoint
x.fixation.start_time = fixation_time(:,1); % timestamp of fixation start
x.fixation.end_time = fixation_time(:,2); % timestamp of fixation end
x.fixation.position = avgpos; % avg. x and y position of fixation
x.blink.start_time = blink_start_time; % timestamp of eye blinks
x.blink.end_time = blink_end_time; % timestamp of eye blinks
% reshape square times to Nx5
x.trial.square_time = reshape(display_square_time,5,numel(x.trial.start_time))'; % time when each square was displayed

% Clean up
clear button_* *_time *pos trial_type % don't clear sync_events - we need it for sync check later
disp('...Success!')


%% GET SQUARE COLOR INFO FROM TRIAL TYPES
disp('Getting square color info...')
target_color = find_first_word(eye_filename,'TargetColor','%s');
if isempty(target_color)
    disp('target color not found - defaulting to Blue!');
    x.target_color = 'B';
else
    x.target_color = target_color(1);
end

middle_color = repmat('-',numel(x.trial.type),5);
% is_right_cross = false(size(x.trial.type));
for i=1:numel(x.trial.type)
    trialstring = num2str(x.trial.type(i));
%     if trialstring(1)=='2'
%         is_right_cross(i) = true;
%     end
    for j=1:5
        if ismember(trialstring(j),'35')
            middle_color(i,j) = 'R';
        elseif ismember(trialstring(j),'16')
            middle_color(i,j) = 'G';
        elseif ismember(trialstring(j),'24')
            middle_color(i,j) = 'B';
        end
    end
end
% x.trial.is_right_cross = is_right_cross;
x.trial.middle_color = middle_color;
x.trial.is_target_color = middle_color==x.target_color;
if strcmp(prefix,'sf3')
    x.trial.is_target_trial = sum(x.trial.is_target_color,2)>2;
else
    x.trial.is_target_trial = sum(x.trial.is_target_color,2)>1;
end

% Clean up
clear middle_color is_right_cross i j trialstring
disp('...Success!')


%% PERFORM DRIFT CORRECTION
figure(100+ascSession);
if drift_correct
    disp('Peforming drift correction...')
    [x,drift,rotation] = perform_drift_correction(x,3,1);
else
    disp('Calculating drift...')
    [~,drift,rotation] = perform_drift_correction(x,3,1);
end
x.trial.drift = drift;
x.trial.rotation = rotation;
clear offset
disp('...Success!')


%% GET OBJECT EVENTS
disp('Getting Saccade Events...')
x.pixel_threshold = pixel_threshold;

target_squares_sofar = cumsum(x.trial.is_target_color,2);
% target_squares_sofar(x.trial.is_right_cross,:) = ...
%     fliplr(cumsum(fliplr(x.trial.is_target_color(x.trial.is_right_cross,:)),2));
x.trial.target_squares_sofar = target_squares_sofar;

% Get trial number for each saccade and fixation
saccade_trialnum = nan(size(x.saccade.start_time));
fixation_trialnum = nan(size(x.fixation.start_time));
for i=1:numel(x.trial.start_time)
    isInTrial = x.saccade.start_time>x.trial.fix_time(i) & x.saccade.start_time>=x.trial.circle_time(i);
    saccade_trialnum(isInTrial) = i;
    isInTrial = x.fixation.start_time>x.trial.fix_time(i) & x.fixation.start_time>=x.trial.circle_time(i);
    fixation_trialnum(isInTrial) = i;
end
x.saccade.trialnum = saccade_trialnum;
x.fixation.trialnum = fixation_trialnum;

% [saccade_trialnum, saccade_squarenum, saccade_class] = classify_squares_saccades(x);
% x.saccade.trialnum = saccade_trialnum;
% x.saccade.squarenum = saccade_squarenum;
% x.saccade.class = saccade_class;
% 
% [fixation_trialnum, fixation_squarenum, fixation_class] = classify_squares_saccades(x,x.fixation.start_time,x.fixation.position);
% x.fixation.trialnum = fixation_trialnum;
% x.fixation.squarenum = fixation_squarenum;
% x.fixation.class = fixation_class;

clear saccade_* fixation_* *squares*
disp('...Success!')


%% GET RESPONSE INFO FROM BUTTON DATA
disp('Getting response info...')

response_time = zeros(size(x.trial.start_time));
reaction_time = zeros(size(x.trial.start_time));
is_target_response = nan(size(x.trial.start_time));
button_trialnum = nan(size(x.button.time));
switch x.experiment       
    case {'SquaresFix_v1pt2' 'SquaresFix3_v1pt0'}
        for i=1:numel(x.trial.start_time)

            iFirstButton = find(x.button.time>x.trial.start_time(i),1);
            if ~isempty(iFirstButton) && i==numel(x.trial.start_time) || x.button.time(iFirstButton)<x.trial.start_time(i+1) % CHANGED FROM V1PT4
                button_trialnum(iFirstButton) = i;
                response_time(i) = x.button.time(iFirstButton);
                reaction_time(i) = response_time(i)-x.trial.circle_time(i); % CHANGED FROM SQUARES V1PT4
                is_target_response(i) = (x.button.number(iFirstButton)==Constants.TARGET_BUTTON);        
            else
                response_time(i) = NaN;
                reaction_time(i) = NaN;
            end    
        end
end
x.button.trialnum = button_trialnum;
x.trial.response_time = response_time;
x.trial.reaction_time = reaction_time;
x.trial.is_target_response = is_target_response;
x.trial.is_correct_response = x.trial.is_target_response==x.trial.is_target_trial;

% Clean up
clear median_iti *_time is_target_response i iFirstButton ct
disp('...Success!')

if ~isempty(eegSession)
    %% GET EEG EVENTS
    disp('Getting EEG events...')
    Constants = GetSquaresConstants;
    x.eeg.system = eegFiletype;
    switch eegFiletype
        case 'Biosemi'
            % Get event sample rate (typically 2048 Hz for Biosemi)
            x.eeg.eventsamplerate = get_biosemi_eventsamplerate(eeg_filename);
            % Get events        
            eeg_sync_events = get_biosemi_events(eeg_filename); % Get raw events
        case {'Sensorium-2008' 'Sensorium-2011' 'Sensorium-2011-1000Hz' 'Sensorium-2011-50Hz'}
            x.eeg.eventsamplerate = 1000;
            % Get events    
            eeg_sync_events = get_sensorium_events(eeg_filename); % Get raw events
        case 'Sensorium-2013'
            x.eeg.eventsamplerate = 1000;
            eeg_sync_data = readEEG_sensorium2013(eeg_filename,1,98); % channel 1 is events
            iEvents = find(diff(eeg_sync_data)~=0)+1;
            eeg_sync_events = [iEvents', eeg_sync_data(iEvents)'];
    end
    if find(eeg_sync_events(:,2)==Constants.BLOCK_START,1)
        eeg_sync_events = eeg_sync_events(find(eeg_sync_events(:,2)==Constants.BLOCK_START,1):end,:); % should start from the 'START_RECORDING' broadcast
    end
    eeg_sync_events = eeg_sync_events(ismember(eeg_sync_events(:,2),[Constants.BLOCK_START, Constants.FIXCROSS_ON, Constants.TRIAL_START, Constants.TRIAL_END, Constants.ENDCIRCLE_OFF, 0]), :);

    [x.sync, eeg_sync_events] = FixSyncEvents(x,eeg_sync_events);
    x.sync.eeg = eeg_sync_events(:,1); % should start from the 'START_RECORDING' broadcast

    disp('...Success!');
    % don't clear sync_events - we need it for sync check later

    %% CONVERT EYELINK EVENTS TO EEG-BASED TIMES

    % disp('Converting eyelink events to eeg times..');
    % 
    % x = EyelinkToEegTimes(x);
    % disp('...Success!');

    %% COMPARE EYE AND EEG SYNC EVENTS
    disp('Checking Parallel Port Events (Eye Tracker vs. EEG)...');

    figure(200+ascSession); clf;
    y.eyelink.sync_events = [x.sync.eyelink, x.sync.events];
    y.eeg.sync_events = eeg_sync_events;
    y.eeg.eventsamplerate = x.eeg.eventsamplerate;
    CheckSyncEvents(y);

    clear y
    disp('Success!');
else
    disp('SKIPPING eeglab events...');
end

%% DISPLAY STATS
disp('Getting Stats...')
disp('---STATS---')
fprintf('%d target trials, %d distractor trials\n',...
    sum(x.trial.is_target_trial),sum(~x.trial.is_target_trial));
% pcomp = mean(~isnan(x.trial.completion_time))*100;
% pcompT = mean(~isnan(x.trial.completion_time(x.trial.is_target_trial)))*100;
% pcompD = mean(~isnan(x.trial.completion_time(~x.trial.is_target_trial)))*100;
% fprintf('%.1f%% trials completed (%.1f%% target, %.1f%% distractor)\n',...
%     pcomp,pcompT,pcompD);
pc = mean(x.trial.is_correct_response)*100; % percent correct
pnr = mean(isnan(x.trial.response_time))*100; % percent no response
pic = 100-pc-pnr; % percent incorrect
fprintf('%.1f%% correct, %.1f%% incorrect, %.1f%% no response\n',pc,pic,pnr);
% fprintf('%d integration saccades, %d completion saccades, %d extra info saccades, %d backward saccades\n',...
%     sum(x.saccade.class==Constants.INTEGRATION), sum(x.saccade.class==Constants.COMPLETION),...
%     sum(x.saccade.class==Constants.EXTRA), sum(x.saccade.class==Constants.BACKWARD));
% fprintf('%d saccades to nothing\n',sum(x.saccade.class==Constants.OTHER));

figure(300+ascSession);
subplot(2,1,1)
hist(sum(x.trial.is_target_color,2),0:5);
xlabel('# target squares')
ylabel('# trials')
subplot(2,1,2)
n = histc(x.trial.reaction_time,0:50:3000);
bar(0:50:3000,n,'histc');
xlim([0 3000])
xlabel('Reaction time (ms)')
ylabel('# trials')
disp('-----------')

% clear pc pnr pic n
disp('...Success!')

%% SAVE
disp('Saving data...')
save(sprintf('%s-%d-%d.mat',prefix,subject,ascSession),'x');
disp('...DONE!')

