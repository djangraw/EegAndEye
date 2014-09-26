function [xAfter yAfter] = ApplyEyeCalibration(xBefore, yBefore, x)

% Convert eye position from Eyelink coordinates to Unity coordinates.
%
% [xAfter yAfter] = ApplyEyeCalibration(xBefore, yBefore, x)
% xyAfter = ApplyEyeCalibration(xyBefore, x)
% xOut = ApplyEyeCalibration(x)
%
% Inputs:
%   - xBefore and yBefore are 1xn vectors of horizontal and vertical eye 
%     positions (e.g., the raw eye position reported by eyelink).
%   - x is a 3DSearch data structure as imported by Import_3DS_data_v2. It
%     must, specifically, have the x.calibration field intact.
%   - xyBefore is a 2-column matrix in which the first column is horizontal
%     eye positions and the second is vertical eye positions before 
%     calibration.
% Outputs:
%   - xAfter and yAfter are the corresponding 1xn vectors of horizontal and
%     vertical eye position that accurately reflect the subject's eye
%     position in screen coordinates during the Unity experiment.
%   - xyAfter is a 2-column matrix in which the first column is horizontal
%     eye positions and the second is vertical eye positions after
%     calibration.
% 
% Created 10/1/10 by DJ.
% Updated 3/28/11 by DJ - allows two-column [x,y] input and output.
% Updated 5/6/11 by DJ - comments.
% Updated 8/8/11 by DJ - added isCalibrated check to prevent calibrating twice
% Updated 11/8/11 by DJ - avoid error for non-raw inputs in isCalibrated
% Updated 5/13/13 by DJ - added x --> xOut format
% Updated 12/5/13 by DJ - added empty input check

% Empty in, empty out
if isempty(xBefore)
    xAfter = [];
    yAfter = [];
    return;
end

% Detect input format
if nargin==2 && nargout==1 % two-column input and output format
    % properly parse inputs
    x = yBefore;
    yBefore = xBefore(:,2);
    xBefore = xBefore(:,1);
end

% x --> xOut format
if nargin==1 && nargout==1
    % get input
    x = xBefore;
    % Adjust saccade and fixation positions    
    x.eyelink.saccade_positions = ApplyEyeCalibration(x.eyelink.saccade_positions,x);
    x.eyelink.saccade_start_positions = ApplyEyeCalibration(x.eyelink.saccade_start_positions,x);
    x.eyelink.fixation_positions = ApplyEyeCalibration(x.eyelink.fixation_positions,x);
    % set output
    xAfter = x;
    return;
end

% Check if this has already been done!
if isCalibrated(xBefore,yBefore,x)
    disp('Eye position has already been calibrated!  Skipping this step...')
    xAfter = xBefore;
    yAfter = yBefore;
else
    

    % Make manual adjustments to calibration!
    if x.subject==2 && ismember(x.session,41:48)
        x.calibration.eye_offset_x = -85;
        x.calibration.eye_gain_x = .96;
    end


    % perform calibration
    xAfter = (xBefore - x.calibration.eye_offset_x) * x.calibration.eye_gain_x;
    yAfter = (yBefore - x.calibration.eye_offset_y) * x.calibration.eye_gain_y;

    if isCalibrated(xAfter,yAfter,x)
        disp('Calibration completed successfully!');
    else
        warning('Calibration does not match x.eyelink.saccade_positions...')
        disp('Calibration complete!')
    end
end


% Detect output format
if nargin==2 && nargout==1 % two-column input and output format
    % concatenate outputs
    xAfter = [xAfter, yAfter];
end

end


% Check whether calibration has been done
function isDone = isCalibrated(xpos,ypos,x)

    if ~isfield(x.eyelink,'saccade_times') || (x.eyelink.saccade_times(1) - x.eyelink.record_time + 1)>length(xpos)
        isDone = false;
        disp('Skipping isCalibrated check...')
    else
        first_saccade_time = x.eyelink.saccade_times(1) - x.eyelink.record_time + 1;
        first_saccade_truth = x.eyelink.saccade_positions(1,:);

        first_saccade_pos = [xpos(first_saccade_time), ypos(first_saccade_time)];
        if first_saccade_pos == first_saccade_truth
            isDone = true;
        else
            isDone = false;
        end
    end

end