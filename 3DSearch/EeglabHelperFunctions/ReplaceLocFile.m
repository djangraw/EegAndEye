function ReplaceLocFile(subject,sessions,newlocfile)

% Replaces the electrode locations on 3DS files with those from a new .loc 
% file and resaves it.
%
% ReplaceLocFile(subject,sessions,newlocfile)
%
% INPUTS:
% -subject is the number of the subject.
% -sessions is either a) a vector of integers indicating the session number
%                        (3DS-<subject>-<sessions>-filtered.set)
%                     b) a string indicating the file suffix 
%                        (3DS-<subject>-<sessions>.set)
%                     c) a vector of cells, each of which contains a string 
%                        as described in b). 
% -newlocfile is the name of the location file, which should be in folder
%  Tools/eeglab/locationfiles. [default: 'jen_sensorium87_79chan.loc']
%
% Created 3/15/11 by DJ.

if nargin<3
    newlocfile = 'jen_sensorium87_79chan.loc';
end

data_dir = [cd '/'];
ALLEEG = [];
newloc_dir = '/Users/dave/Documents/Tools/eeglab/locationfiles/';
if ischar(sessions) % for a string indicating the file suffix, sessions must have length 1 to work in the loop below
    sessions = {sessions}; % so place the string in a cell.
end

for i=1:numel(sessions)
    % load file
    if isnumeric(sessions)
        filename = sprintf('3DS-%d-%d-filtered.set',subject,sessions(i));
    elseif iscell(sessions)
        filename = sprintf('3DS-%d-%s.set',subject,sessions{i});
    end
    disp(['---Replacing loc file for ' filename '...'])
    EEG = pop_loadset('filename',filename,'filepath',data_dir);
    % load new loc file
    EEG=pop_chanedit(EEG, 'load',{[newloc_dir newlocfile] 'filetype' 'autodetect'});
    EEG = eeg_checkset( EEG );
    % save
    EEG = pop_saveset( EEG, 'savemode','resave');
    disp('   ...Done.')
end