function [ALLEEG, EEG, CURRENTSET] = InterpolateElectrodes(input_filename, output_filename, chanList)

% Loads a dataset, interpolates the specified electrodes, saves the result.
%
% [ALLEEG, EEG, CURRENTSET] = InterpolateElectrodes(input_filename,
% output_filename, chanList)
%
% INPUTS:
% -input_filename is the filename of an existing eeglab file in the current
% directory - the one you want to remove electrodes from.
% -output_filename is the filename where you want to save the resulting
% file.
% -chanList is a vector of channel numbers you wish to interpolate, or a
% vector of cells containing the names of the channels you wish to interp.
% 
% OUTPUTS:
% -ALLEEG, EEG, and CURRENTSET are the standard eeglab variables of the
% same names.
%
% Created 8/8/12 by DJ based on RemoveElectrodes.m.
% Updated 3/18/13 by DJ - automatically adds '.set' if it's not in filename
% Updated 2/19/15 by DJ - switch to EEGLAB's spherical interpolation

%% Create inputs and set up
% Check filename
if length(input_filename)<4 || ~strcmp(input_filename(end-3:end),'.set'); % if the data file doesn't end in .set...
    input_filename = strcat(input_filename,'.set'); % add on .set
end
% Data info
data_dir = [cd '/'];
% data_dir = '/Users/dave/Documents/Data/3DSearch/';

fprintf('---Loading file %s, interpolating %d channels, and saving to %s---\n',input_filename, numel(chanList),output_filename);
ALLEEG = [];

%% Load dataset
% disp('Loading...')
EEG = pop_loadset('filename',input_filename,'filepath',data_dir);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

%% Get channel numbers (added 2/19/15)
% Convert input to vector of channel indices
if isstruct(chanList)
    [~, chanNumVec] = ismember({chanList.labels},{EEG.chanlocs.labels});
elseif iscell(chanList)
    [~, chanNumVec] = ismember(chanList,{EEG.chanlocs.labels});
elseif isnumeric(chanList)
    chanNumVec = chanList;
end

%% Remove channels
% disp('Removing channels...')
% EEG = InterpChan(EEG,chanList); % remove electrodes
EEG = eeg_interp(EEG,chanNumVec); % use EEGLAB function with spherical interpolation
EEG = pop_editset(EEG, 'setname', sprintf('%s - %dChansInterpolated',EEG.setname,numel(chanList))); % rename dataset
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 

%% Save results
% disp('Saving...')
EEG = pop_saveset( EEG, 'filename',output_filename, 'filepath',data_dir);
disp('Data has been saved!');
