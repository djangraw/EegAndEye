function EEG = rereferenceData(EEG,refChans)

% Re-references the data currently selected in EEGLAB to the channels
% specified in the variable refChans. (default: no re-referencing).
%
% EEG = rereferenceData(EEG,refChans)
%
% INPUTS:
% -EEG is an eeglab data structure.
% -refChans can be a vector of channel numbers, a vector of cells
% containing strings indicating channels, or one of the following string 
% code-words: 'Average', 'Cz','TP7/TP8','Non-frontal'.
%
% OUTPUT:
% -EEG is the inputted eeglab data structure, rereferenced as specified.
%
% Created 11/2/10 by DJ.
% Updated 11/4/10 by DJ - allowed different inputs for chans
% Updated 11/8/10 by DJ - changed dataset naming convention
% Updated 2/21/11 by DJ - allowed refChans to be an input (sort of), and
% will now operate if no re-referencing is desired.
% Updated 2/28/11 by DJ - made a function
% Updated 6/2/11 by DJ - commented out ALLEEG line (ALLEEG not declared)

%% Set up
if nargin<2 || isequal(refChans,'None');
    disp('Not re-referencing data...')
    return;
end

switch refChans
    case 'Average'
        refChans = {''};
    case 'Cz'
        refChans = {'Cz'};
    case 'TP7/TP8'
        refChans = {'TP7', 'TP8'};
    case 'Non-frontal' % all except the front-most three rows
        refChans = {'FT7', 'FC5', 'FC3', 'FC1', 'C1', 'C3', 'C5', 'T7', 'TP7', ...
            'CP5', 'CP3', 'CP1', 'P1', 'P3', 'P5', 'P7', 'P9', 'PO7', 'PO3', ...
            'O1', 'Iz', 'Oz', 'POz', 'Pz', 'CPz', 'FT8', 'FC6', 'FC4', 'FC2', ...
            'FCz', 'Cz', 'C2', 'C4', 'C6', 'T8', 'TP8', 'CP6', 'CP4', 'CP2', ...
            'P2', 'P4', 'P6', 'P8', 'P10', 'PO8', 'PO4', 'O2' };
end

% convert from string to numbers, if necessary
chanLabels = {EEG.chanlocs(:).labels}; % the names of all the channels
if iscell(refChans)
    refChans = find(ismember(chanLabels,refChans));
end

%% Re-reference to average of specified channels
fprintf('Re-referencing to channel(s) [%s]...\n',num2str(refChans))
EEG = pop_reref( EEG, refChans,'keepref','on');
EEG = pop_editset(EEG, 'setname', [EEG.setname '-rereferenced']); % rename dataset
EEG = eeg_checkset( EEG );
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% Clean up
% clear refChans chanLabels