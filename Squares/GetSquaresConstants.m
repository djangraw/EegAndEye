function Constants = GetSquaresConstants(when)

% Constants = GetSquaresConstants(when)
%
% Created 10/19/11 by DJ.
% Updated 10/26/11 by DJ - changed button constants
% Updated 10/27/11 by DJ - added base and blink constants
% Updated 10/31/11 by DJ - changed SACCADE related constants
% Updated 11/28/11 by DJ - added lock constants, eventnames field
% Updated 8/1/12 by DJ - changed blinks to 'BS', 'BE'
% Updated 3/5/13 by DJ - added SQUARE_ON

if nargin<1
    when = now;
end

% TTL Pulse Constants
Constants.BLOCK_START = 30;
Constants.FIXCROSS_ON = 40;
Constants.TRIAL_START = 50;
Constants.SQUARE_ON = 55; % for squaresfix experiment
Constants.TRIAL_END = 60;
Constants.ENDCIRCLE_OFF = 70;
Constants.ISTARGET = 1; % add this for target trials in EEGLAB files

% Button Constants
Constants.BUTTON_BASE = 100; % events in EEGLAB files will start here
Constants.TARGET_BUTTON = 3; % eyelink button number for target response
Constants.NONTARGET_BUTTON = 1; % eyelink button number for nontarget response

% Square Position Constants
Constants.SQUARE_X = [212, 362, 512, 662, 812];
Constants.SQUARE_Y = repmat(384,1,5);
Constants.LEFTCROSS_X = 62;
Constants.LEFTCROSS_Y = 384;
Constants.RIGHTCROSS_X = 962;
Constants.RIGHTCROSS_Y = 384;

% Saccade type constants
Constants.SACCADESTART_BASE = 120; % event codes in EEGLAB will start here
Constants.SACCADEEND_BASE = 150;
Constants.DISTRACTOR = 0; % Saccade to a distractor square
Constants.INTEGRATION = 1; % Saccade to first target square
Constants.COMPLETION = 2; % Saccade to second target square
Constants.EXTRA = 3; % Saccade to third/fourth/fifth target square
Constants.WITHIN = 4; % Saccade made within a single square
Constants.BACKWARD = 5; % Saccade from a square to the one before it
Constants.FIXCROSS = 6; % Saccade to fixation cross
Constants.ENDCIRCLE = 7; % Saccade to end circle
Constants.OTHER = 19; % Saccade to no known object

% Fixation direction constants
Constants.FIXSTART_BASE = 180;
Constants.FIXEND_BASE = 190;
Constants.LEFTSIDE = 0;
Constants.RIGHTSIDE = 1;
Constants.OTHERSIDE = 2;

% Other EEGLAB event constants
Constants.BLINKSTART = 200;
Constants.BLINKEND = 201;

% Epoch Locking Events
Constants.DISTRACTOR_LOCK = 250;
Constants.TARGET_LOCK = 260;

% Event Names
Constants.EVENTNAMES{Constants.BLOCK_START} = 'BlockStart';
Constants.EVENTNAMES{Constants.FIXCROSS_ON+0} = 'FixOn-D';
Constants.EVENTNAMES{Constants.FIXCROSS_ON+1} = 'FixOn-T';
Constants.EVENTNAMES{Constants.TRIAL_START+0} = 'TrialStart-D';
Constants.EVENTNAMES{Constants.TRIAL_START+1} = 'TrialStart-T';
Constants.EVENTNAMES{Constants.SQUARE_ON} = 'SquareOn';
Constants.EVENTNAMES{Constants.TRIAL_END+0} = 'TrialEnd-D';
Constants.EVENTNAMES{Constants.TRIAL_END+1} = 'TrialEnd-T';
Constants.EVENTNAMES{Constants.BUTTON_BASE+0} = 'Button-D';
Constants.EVENTNAMES{Constants.BUTTON_BASE+1} = 'Button-T';
Constants.EVENTNAMES{Constants.BUTTON_BASE+2} = 'Button-X'; % other
Constants.EVENTNAMES{Constants.ENDCIRCLE_OFF+Constants.LEFTSIDE} = 'FixOff-L';
Constants.EVENTNAMES{Constants.ENDCIRCLE_OFF+Constants.RIGHTSIDE} = 'FixOff-R';
Constants.EVENTNAMES{Constants.SACCADESTART_BASE} = 'SS';%'SacStart';
Constants.EVENTNAMES{Constants.SACCADEEND_BASE} = 'SE';%'SacEnd';
Constants.EVENTNAMES{Constants.FIXSTART_BASE+Constants.LEFTSIDE} = 'FSL';%'FixStart-L';
Constants.EVENTNAMES{Constants.FIXEND_BASE+Constants.LEFTSIDE} = 'FEL';%'FixEnd-L';
Constants.EVENTNAMES{Constants.FIXSTART_BASE+Constants.RIGHTSIDE} = 'FSR';%'FixStart-R';
Constants.EVENTNAMES{Constants.FIXEND_BASE+Constants.RIGHTSIDE} = 'FER';%'FixEnd-R';
Constants.EVENTNAMES{Constants.BLINKSTART} = 'BS';
Constants.EVENTNAMES{Constants.BLINKEND} = 'BE';
Constants.EVENTNAMES{Constants.DISTRACTOR_LOCK} = 'TimeLock-D';
Constants.EVENTNAMES{Constants.TARGET_LOCK} = 'TimeLock-T';
Constants.EVENTNAMES{Constants.SQUARE_ON} = 'Square';
