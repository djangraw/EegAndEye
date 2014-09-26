function EEG = CombineEeglabSessions(subject,sessions,inputSuffix,outputSuffix,experimentType)

% Combines all the eeglab files into one long file.  This should be done 
% after all filtering is done (ImportToEeglab) and events are loaded 
% (AddEeglabEvents), but before epoching.
%
% EEG = CombineEeglabSessions(subject,sessions,inputSuffix,outputSuffix)
%
% - EEGLAB should already be started.
%
% INPUTS:
% -subject and session are numbers that tell us which file to save
% -inputSuffix is a string specifying the end of the single-file filenames.
% Full filename is '3DS-<subject>-<session><inputSuffix>.set' [-filtered]
% -outputSuffix is a string specifying the end of the combined filename.  
% Full filename is '3DS-<subject><inputSuffix>.set' [-all-filtered]
% -experimentType is a string specifying the type of experiment ('3DSearch,
% 'Squares','TargetSwitch') or the prefix to the filename (before 
% '-<subject>').
%
% OUTPUT:
% -EEG is the combined file's eeglab struct.
%
% - TO DO: check on end-of-session problems with epoching?
% 
% Created 8/2/10 by DJ.
% Updtaed 8/20/10 by DJ - comments
% Updated 2/23/11 by DJ - made a function. 
% Updated 3/1/11 by DJ - comments
% Updated 10/27/11 by DJ - added experimentType input to work with Squares
%  experiments
% Updated 10/31/11 by DJ - keep record of merged files in 
%  EEG.etc.mergedfiles
% Updated 5/3/12 by DJ - added TargetSwitch support, experimentType->prefix
% Updated 7/3/12 by DJ - added check that there are >1 sessions
% Updated 3/17/13 by DJ - added SquaresFix3 compatibility.

%% CHECK INPUTS AND SET UP
if nargin<3 
    inputSuffix = '-filtered';
end
if nargin<4 
    outputSuffix = '-all-filtered';
end
if nargin<5 || strcmp(experimentType,'');
    experimentType = '3DSearch';
end

% Set up
switch experimentType
    case '3DSearch'
        prefix = '3DS';        
    case 'Squares'
        prefix = 'sq';
    case 'TargetSwitch'
        prefix = 'TS';
    case 'SquaresFix'
        prefix = 'sf';
    case 'SquaresFix3'
        prefix = 'sf3';
    otherwise
        prefix = experimentType;
end
data_dir = [cd '/'];
% data_dir = '/Users/dave/Documents/Data/3DSearch/';
final_filename = sprintf('%s-%d%s',prefix,subject,outputSuffix);

% Clear datasets from EEGLAB
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
eeglab redraw;

fprintf('---Combining data for subject %d (%d sessions) and saving to %s---\n',subject,numel(sessions),final_filename);
tic

%% LOAD DATASETS
for i=1:numel(sessions)    
    filename{i} = sprintf('%s-%d-%d%s.set',prefix,subject,sessions(i),inputSuffix);
    EEG = pop_loadset('filename',filename{i},'filepath',data_dir);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
end

%% MERGE DATASETS
if numel(sessions)>1
    EEG = pop_mergeset( ALLEEG, [1:numel(sessions)], 0);
end
EEG.etc.mergedfiles = filename;
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',sprintf('%s-%d-all',prefix,subject),'gui','off'); 
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% SAVE DATASET
EEG = pop_saveset( EEG, 'filename',final_filename,'filepath',data_dir);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% CLEAN UP
toc % Display elapsed time
% eeglab redraw
% clear data_dir final_filename i