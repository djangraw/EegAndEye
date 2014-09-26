function RemoveEogArtifacts(input_filename,output_filename,SACCADE_EXPAND,OUTSIDE_WINDOW)

% Given an unepoched EEGLAB data file with blink and saccade events 
% included, find the maximum difference component within all data inside saccades,
% project it out of that saccade data, and save the result.
% 
% EEG = RemoveEogArtifacts(input_filename,output_filename,SACCADE_EXPAND,OUTSIDE_WINDOW)
%
% INPUTS:
% -input_filename is the filename of an existing eeglab file in the current
% directory - the one you want to remove eog from.
% -output_filename is the filename where you want to save the resulting
% file.
% - SACCADE_EXPAND is a scalar indicating how far eog artifacts extend 
% beyond the boundaries defined by eyelink events.
% - OUTSIDE_WINDOW is a scalar indicating how much data on either side of
% each saccade should be averaged during mean subtraction.  if this = inf,
% the component is removed from the entire dataset.
%
% Created 10/28/11 by DJ based on RemoveBlinkArtifacts.

%% Handle Inputs
if nargin<3
    SACCADE_EXPAND = 0; % time, in ms, beyond the saccade start/end events surrounding a blink, that you want to be considered in the blink.
end
if nargin<4
    OUTSIDE_WINDOW = 25; % time, in ms, on either side of the blink that we want to average to do mean subtraction
end

%% Create inputs and set up
% Data info
data_dir = [cd '/'];
% data_dir = '/Users/dave/Documents/Data/3DSearch/';

fprintf('---Loading file %s, removing eog artifacts, and saving to %s---\n',input_filename, output_filename);

%% Load dataset
disp('Loading...')
EEG = pop_loadset('filename',input_filename,'filepath',data_dir);

%% Remove BlinkData
EEG = eogartifacts(EEG,SACCADE_EXPAND,OUTSIDE_WINDOW);

%% Save results
disp('Saving...')
EEG = pop_saveset( EEG, 'filename',output_filename, 'filepath',data_dir);
disp('Data has been saved!');
