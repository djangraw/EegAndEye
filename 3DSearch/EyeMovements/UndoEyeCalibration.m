function [xAfter, yAfter] = UndoEyeCalibration(xBefore, yBefore, x)

% The opposite of ApplyEyeCalibration.
% 
% [xAfter, yAfter] = UndoEyeCalibration(xBefore, yBefore, x)
% xyAfter = UndoEyeCalibration(xyBefore, x)
% xOut = UndoEyeCalibration(x)
%
% Inputs:
%   - xBefore and yBefore are 1xn vectors of horizontal and vertical eye 
%     positions (after calibration).
%   - x is a 3DSearch data structure as imported by Import_3DS_data_v2. It
%     must, specifically, have the x.calibration field intact.
%   - xyBefore is a 2-column matrix in which the first column is horizontal
%     eye positions and the second is vertical eye positions (after 
%     calibration).
% Outputs:
%   - xAfter and yAfter are the corresponding 1xn vectors of horizontal and
%     vertical eye position (with calibration removed, i.e., the raw eye
%     position reported by eyelink).
%   - xyAfter is a 2-column matrix in which the first column is horizontal
%     eye positions and the second is vertical eye positions (with
%     calibration removed).
%
% Created 5/13/13 by DJ based on ApplyEyeCalibration.

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
    x.eyelink.saccade_positions = UndoEyeCalibration(x.eyelink.saccade_positions,x);
    x.eyelink.saccade_start_positions = UndoEyeCalibration(x.eyelink.saccade_start_positions,x);
    x.eyelink.fixation_positions = UndoEyeCalibration(x.eyelink.fixation_positions,x);
    % Reset calibration struct
    x.calibration.eye_offset_x = 0;
    x.calibration.eye_offset_y = 0;
    x.calibration.eye_gain_x = 1;
    x.calibration.eye_gain_y = 1;
    % set output
    xAfter = x;
    return;
end


% perform calibration
xAfter = (xBefore / x.calibration.eye_gain_x) + x.calibration.eye_offset_x;
yAfter = (yBefore / x.calibration.eye_gain_y) + x.calibration.eye_offset_y;

% Detect output format
if nargin==2 && nargout==1 % two-column input and output format
    % concatenate outputs
    xAfter = [xAfter, yAfter];
end
