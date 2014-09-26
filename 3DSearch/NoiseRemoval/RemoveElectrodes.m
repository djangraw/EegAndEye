function [ALLEEG, EEG, CURRENTSET] = RemoveElectrodes(input_filename, output_filename, chans_toremove)

% RemoveElectrodes
%
% Gets rid of frontal electrodes (or any other electrodes you specify
% the code), then saves the result as a new file.
%
% INPUTS:
% -input_filename is the filename of an existing eeglab file in the current
% directory - the one you want to remove electrodes from.
% -output_filename is the filename where you want to save the resulting
% file.
% -chans_toremove is a vector of channel numbers you wish to remove, or a
% vector of cells containing the names of the channels you wish to remove.
% 
% OUTPUTS:
% -ALLEEG, EEG, and CURRENTSET are the standard eeglab variables of the
% same names.
%
% Created 10/26/10 by DJ. (called RemoveFrontalElectrodes)
% Updated 11/2/10 by DJ - removed re-referencing.
% Updated 11/29/10 by DJ - switched to RemoveElectrodes, made into a fn.
% Updated 11/30/10 by DJ - added Convert to numbers section 
% Updated 3/18/13 by DJ - automatically adds '.set' if it's not in filename
% Updated 2/25/14 by DJ - document which chans were removed

%% Create inputs and set up
% Check filename
if length(input_filename)<4 || ~strcmp(input_filename(end-3:end),'.set'); % if the data file doesn't end in .set...
    input_filename = strcat(input_filename,'.set'); % add on .set
end
% Data info
data_dir = [cd '/'];
% data_dir = '/Users/dave/Documents/Data/3DSearch/';

fprintf('---Loading file %s, removing %d channels, and saving to %s---\n',input_filename, numel(chans_toremove),output_filename);
ALLEEG = [];

%% Load dataset
% disp('Loading...')
EEG = pop_loadset('filename',input_filename,'filepath',data_dir);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

%% Convert to numbers, if necessary
chanLabels = {EEG.chanlocs(:).labels}; % the names of all the channels
if iscell(chans_toremove)
    chans_toremove = find(ismember(chanLabels,chans_toremove));
end

% clear chanLabels

%% Remove channels
% disp('Removing channels...')
EEG = pop_select( EEG,'nochannel',chans_toremove); % remove electrodes
EEG = pop_editset(EEG, 'setname', sprintf('%s - %dChannelsRemoved',EEG.setname,numel(chans_toremove))); % rename dataset
EEG.etc.chansRemoved = chanLabels(chans_toremove); % document which chans were removed
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 

%% Save results
% disp('Saving...')
EEG = pop_saveset( EEG, 'filename',output_filename, 'filepath',data_dir);
disp('Data has been saved!');
