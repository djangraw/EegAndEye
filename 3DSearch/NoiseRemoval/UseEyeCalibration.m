function EEG = UseEyeCalibration(search_filename,eye_filename,final_filename,refChans)

% Loads an eyecalibration data file produced by ImportEyeCalibration.m
% and a regular data file (usually a 3DSearch data file produced by 
% CombineEeglabSessions), and employs the pop_eyesubtract() function. 
%
% EEG = UseEyeCalibration(subject,session,version,refChans)
%
% - This should find the eye blink/motion artifact components in the 
% eyecalibration file and remove them from the regular data file.
% - The 'clean' data file produced will be saved in the same directory that
% the original files were in.
%
% INPUTS:
% - search_filename is an eeglab file containing eye artifact in the 
% current path. (example: '3DS-<subject>-all-filtered.set')
% - eye_filename is an eeglab file containing eyecalibration data in the 
% current path (example:'eyecalibration-<subject>-<session>.set')
% - final_filename is the eeglab file to which we want to save our output
% (example: '3DS-<subject>-all-eyeremoved.set')
%
% OUTPUT:
% - EEG is the eeglab struct of the new, 'clean' dataset.
%
% Created 8/20/10 by DJ.
% Updated 2/21/11 by DJ - adjusted to new rereferenceData input, added
% 'version' input
% Updated 3/2/11 by DJ - made function, moved reref to this function from 
% ImportEyeCalibration().
% TO DO: test, update ImportAll/gui

%% CHECK INPUTS AND SET UP
if nargin<4
    refChans = 'None';
end

data_dir = [cd '/'];
% data_dir = '/Users/dave/Documents/Data/3DSearch/';

% Clear datasets from EEGLAB
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
eeglab redraw;

fprintf('---Applying eyecalibration from %s to %s and saving to %s---\n',eye_filename,search_filename,final_filename);
tic

%% Load the datasets
EEG = pop_loadset('filename',eye_filename,'filepath',data_dir); % eyecalibration
EEG = rereferenceData(EEG,refChans); % Re-reference data as requested
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

EEG = pop_loadset('filename',search_filename,'filepath',data_dir); % 3DSearch
EEG = rereferenceData(EEG,refChans); % Re-reference data as requested
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

%% Apply the calibration
ALLEEG = pop_eyesubtract(ALLEEG, 1, 2, 50, 150, 175, 200, 225, 0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'retrieve',3,'study',0); 

%% Save results
EEG = pop_saveset( EEG, 'filename',final_filename,'filepath',data_dir);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% Clean up
toc % Display elapsed time
eeglab redraw
clear data_dir eye_filename search_filename final_filename

