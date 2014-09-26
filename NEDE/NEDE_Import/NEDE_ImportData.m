function x = NEDE_ImportData(experiment,subject,session,pixelThresholds)

% Loads a NEDE experiment file, compiles the data into a struct, and saves.
%
% x = Import_NEDE_Data(experiment,subject,ascSession,pixelThresholds)
%
% Saves the result in the current directory as 
% '[experiment]-[subject]-[session].mat'.
%
% INPUTS:
% - experiment is a string
% - subject is a scalar
% - session is a scalar
% - pixelThresholds is a vector
%
% OUTPUTS:
% - x is the struct with all the info in it.
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
% Updated 1/8/13 by DJ - v4: removed event_filename info - everything's in the eyelink file! 
% Updated 12/6/13 by DJ - adapted into NEDE version
% Updated 2/19/14 by DJ - comments.
% Updated 9/10/14 by DJ - specify event types to log, added recording start
% Updated 9/25/14 by DJ - fixed bug where pre-trial port messages were ignored


%% SET UP
% Declare defaults
if nargin<4 || isempty(pixelThresholds)
    pixelThresholds = 1000; % include all saccades
end
% Handle multiple sessions
if numel(session)>1
    for i=1:numel(session)
        x(i) = NEDE_ImportData(experiment,subject,session(i),pixelThresholds);
    end
    return;
end


%% GET FILENAME
tic
fprintf('---Importing Data For Subject %d, Session %d---\n',subject,session);
eye_filename = sprintf('%s-%d-%d.asc',experiment,subject,session); % eyelink .edf file converted to text using Visual EDF2ASC.

    
%% READ EYELINK LOG PARAMETERS
disp('Getting Session Parameters...')
x.params = NEDE_ParseParams(eye_filename);
disp('...Success!')

%% READ EVENTS FILE
% This is used for events like objects entering/exiting the scene, 
% saccades to objects, etc.
disp('Getting Object Info...')
objects = NEDE_ParseObjects(eye_filename);
x.objects = objects{1}; % for single-trial experiments
disp('...Success!')

%% READ EYE TRACKER FILE
% This is used for events like objects entering/exiting the scene, 
% saccades, and for checking the eyelink and eeg times against each other.
disp('Getting Session Events...')

% x.events = NEDE_ParseEvents(eye_filename,[],'','END TRIAL');
types = {'saccade','fixation','blink','fixupdate','button','trial','port','camera','visible'};
x.events = NEDE_ParseEvents(eye_filename,types,'','END TRIAL');

% Add messages only from after START TRIAL flag (to avoid counting parameter lines)
foo = NEDE_ParseEvents(eye_filename,{'message'},'START TRIAL','END TRIAL');
x.events.message = foo.message;

% Add recording start time to analysis
foo = NEDE_ParseEvents(eye_filename,{'START'},'','PRESCALER');
x.events.recording = foo.recording;

% Convert eye position from eyelink coordinates to Unity coordinates
fprintf('Calibrating saccade start positions...')
x.events.saccade.position_start = NEDE_ApplyEyeCalibration(x.events.saccade.position_start, x.params.eyelink);
fprintf('Calibrating saccade end positions...')
x.events.saccade.position_end = NEDE_ApplyEyeCalibration(x.events.saccade.position_end, x.params.eyelink);
fprintf('Calibrating fixation positions...')
x.events.fixation.position = NEDE_ApplyEyeCalibration(x.events.fixation.position, x.params.eyelink);
% NOTE: Fixupdates have already had the calibration applied, so we don't
% need to do it here.

% Clean up
disp('...Success!')

%% GET OBJECT EVENTS
disp('Getting Saccade Events...')
if ~isempty(x.events.saccade.position_end)
    [objectSeen,isFirstToObject,isLastToObject] = NEDE_ClassifySaccades(x.events.saccade,x.events.visible,pixelThresholds);
    x.events.saccade.object_seen = objectSeen;
    x.events.saccade.isFirstToObject = isFirstToObject;
    x.events.saccade.isLastToObject = isLastToObject;
else
    [x.events.saccade.object_seen, x.events.saccade.isFirstToObject, x.events.saccade.isLastToObject] = deal([]);
end
% Clean up
disp('...Success!')

%% SAVE STRUCT
disp('Saving result in current directory...')
save(sprintf('%s-%d-%d.mat',experiment,subject,session),'x');
disp('...Success!')
toc % Display elapsed time