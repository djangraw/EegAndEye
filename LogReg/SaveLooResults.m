function SaveLooResults(EEG1,EEG2, time, Az, trainingwindowlength, trainingwindowinterval, datetime, reftype, bootstrap)

% Saves the results of a leave-one-out cross-validation in the file
% 'LOO.mat'.
%
% SaveLooResults(EEG1, EEG2, time, Az, trainingwindowlength, trainingwindowinterval, datetime)
%
% INPUTS:
% -EEG is the EEGLAB file
% -time is a vector of times at the center of the training windows
% -Az is a vector of the leave-one-out classification Az values for each
% window.
% -trainingwindowlength is the length of each training window (in samples).
% -trainingwindowinterval is the distance between each training window
% center. [unique(diff(time))]
% -datetime is the date and time at which the results were saved. [datestr(now)]
% -reftype is a cell array containing strings specifying which channels were used as a reference. [EEG.ref]
% -bootstrap is a binary variable saying whether the result used
% bootstrapping.
%
% Created 11/4/10 by DJ.
% Updated 11/8/10 by DJ - added EEG2 input, setnames
% Updated 1/14/11 by DJ - added bootstrap
% Updated 2/23/11 by DJ - added displays

%% SET UP
% handle inputs
if nargin<9 || isempty(bootstrap)
    bootstrap = 0;
end
if nargin<8 || isempty(reftype)
    reftype = EEG1.ref;
end
if nargin<7 || isempty(datetime)
    datetime = datestr(now);
end
if nargin<6 || isempty(trainingwindowinterval)
    trainingwindowinterval = unique(diff(time));
end

% load loo results
if bootstrap
    looLog = which('LOO_bootstrap.mat'); % finds 'LOO_bootstrap.mat' in the current path
    disp('Loading log LOO_bootstrap.mat...')
else
    looLog = which('LOO.mat'); % finds 'LOO.mat' in the current path
    disp('Loading log LOO.mat...')
end

load(looLog); % loads the variable 'LOO'

%% ADD INFO
LOO = [LOO, struct('setname1','','setname2','','reference','','time',[],'Az',[],'trainingwindowlength',[],'trainingwindowinterval',[],'datetime','')];
LOO(end).setname1 = EEG1.setname;
LOO(end).setname2 = EEG2.setname;
LOO(end).reference = reftype;
LOO(end).time = time;
LOO(end).Az = Az;
LOO(end).trainingwindowlength = trainingwindowlength;
LOO(end).trainingwindowinterval = trainingwindowinterval;
LOO(end).datetime = datetime;

%% SAVE
disp('Saving LOO log...')
save(looLog,'LOO'); % saves the variable
disp('Success!')