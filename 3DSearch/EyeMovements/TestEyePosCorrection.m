function calibration = TestEyePosCorrection(subject,sessions,pixelThresholds,timeLimits)

% calibration = TestEyePosCorrection(subject,sessions,pixelThresholds,timeLimits)
%
% INPUTS:
% -subject is a scalar indicating the number of the 3DSearch subject.
% -sessions is an N-element vector indicating the session numbers you're 
% interested in.
% -pixelThresholds is a scalar indicating the distance from the object that
% will be considered a saccade to the object.
% -timeLimits is a 2-element vector indicating the time range (in ms
% relative to the time when an object appeared) in which a saccade could be
% to the object.
%
% OUTPUTS:
% -calibration is an N-element vector of structs that can be substituted in
% for the x.calibration struct in a 3ds data file.
%
% Created 5/9/13 by DJ.

% Handle inputs
if nargin<3
    pixelThresholds = 200;
end
if nargin<4
    timeLimits = [0 Inf];
end

% declare constants
n = numel(sessions);
CTR_X = 1024/2; % x pos of center of screen
CTR_Y = 768/2; % y pos of center of screen

% Get median positions
[~,medianpos] = EyePositionHeatMap(subject,sessions,'3DS');

% Main loop
for i=1:n
    % Create new calibration
    foo = load(sprintf('3DS-%d-%d.mat',subject,sessions(i)));
    x = UndoEyeCalibration(foo.x);
    calibration(i).eye_offset_x = medianpos(i,1)-CTR_X;
    calibration(i).eye_offset_y = medianpos(i,2)-CTR_Y;
    calibration(i).eye_gain_x = 1;
    calibration(i).eye_gain_y = 1;
    x.calibration = calibration(i);
    fprintf('subj %d, session %d: x = %.1f, y = %.1f\n',subject,sessions(i),x.calibration.eye_offset_x,x.calibration.eye_offset_y);
    
    % apply eye calibration
    x.eyelink.saccade_positions = ApplyEyeCalibration(x.eyelink.saccade_positions, x);
    x.eyelink.saccade_start_positions = ApplyEyeCalibration(x.eyelink.saccade_start_positions, x);
    x.eyelink.fixation_positions = ApplyEyeCalibration(x.eyelink.fixation_positions, x);
    % re-calculate saccade events
    x.eyelink.saccade_events = classify_saccades(x.eyelink.saccade_times,x.eyelink.saccade_positions,x.eyelink.object_limits,pixelThresholds,timeLimits);
    x = EyelinkToEegTimes(x); % to reset the x.eeg.saccade events used in MakeEyeMovie
    
    % Make eye movie
    foo = load(sprintf('3DS-%d-%d-eyepos.mat',subject,sessions(i)));
    MakeEyeMovie(foo.eyepos,foo.pupilsize,x);    
end