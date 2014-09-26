function Numbers = NEDE_GetNumbers

% Get the constants specific to the 3DSearch task. 
%
% Numbers = NEDE_GetNumbers
% GetNumbers; (this is equivalent to Numbers = NEDE_GetNumbers;)
%
% INPUTS:
%
% OUTPUTS:
% -Numbers is the struct of NEDE constants.
%
% Should be similar to the constants section of Numbers.js from the 
% current experiment (adjusting variable declarations for MATLAB).
%
% Function GetNumberMeaning.m can be used to display the meaning of an
% event code you encounter.  When you update this function, you should
% consider updating that one too.
%
% Created 6/10 by DJ
% Updated 7/29/10 by DJ (added Target/Distractor codes)
% Updated 8/19/10 by DJ (added threshold values)
% Updated 10/18/10 by DJ (added BLINK).
% Updated 2/24/11 by DJ (added SEES, SIM_BUTTON).
% Updated 4/18/11 by DJ (multiplied all anchoring constants by 10 to
% accommodate new levels with >50 objects).
% Updated 6/2/11 by DJ - added _START and _END codes for saccades,
%   fixations, and smooth pursuit
% Updated 8/4/11 by DJ - added optional input 'when' (for old files)
% Updated 12/5/13 by DJ - simplified for NEDE release.


%CONSTANTS FOR PARALLEL PORT MESSAGES (MUST BE 0-255).
Numbers.START_RECORDING = 200;
Numbers.END_RECORDING = 201;
Numbers.SYNC = 211;

% Trial Type values:
Numbers.STATIONARY = 0;
Numbers.FOLLOW = 1;

% Button values:
Numbers.GAMEPAD_Y = 1;
Numbers.GAMEPAD_X = 2;
Numbers.GAMEPAD_B = 3;
Numbers.GAMEPAD_A = 4;
Numbers.GAMEPAD_ARROWS = 5;
Numbers.GAMEPAD_L = 6;
Numbers.GAMEPAD_R = 7;

% Button assignments:
Numbers.GAMEPAD_TARGET_BUTTON = Numbers.GAMEPAD_B; % number of the target (A) button on the gamepad (?)
Numbers.GAMEPAD_BRAKE_BUTTON = Numbers.GAMEPAD_A; % number of the brake (B) button on the gamepad

% ANALYSIS CONSTANTS
Numbers.MAX_LAG_BETWEEN_FRAMES = 100; % max time (ms) something can be offscreen before it's considered a real disappearance

    
% Send results to program that called this one
if nargout==0
    assignin('caller','Numbers',Numbers);
end