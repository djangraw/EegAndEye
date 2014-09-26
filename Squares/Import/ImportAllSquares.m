% ImportAllSquares script

% Created 10/27/11 by DJ.
% Updated 10/28/11 by DJ - added post-import functions

%% SET UP
subject = 3;
ascSessions = 1:6;
eegSessions = ascSessions-1;

%% IMPORT DATA
for i=1:numel(ascSessions)
    import_squares_data(subject,ascSessions(i),eegSessions(i),'Sensorium-2011',50);
    ImportToEeglab(subject,ascSessions(i),eegSessions(i),'Sensorium-2011','Squares');
    AddEeglabEvents(subject,ascSessions(i),'-filtered','Squares');
end
CombineEeglabSessions(subject,ascSessions,'-filtered','-all-filtered','Squares');

%% REMOVE BAD ELECTRODES

inputsuffix = '-all-filtered';
outputsuffix = '-all-noduds';
duds = {'NZ', 'F9', 'F10'};
RemoveElectrodes(sprintf('sq-%d%s.set',subject,inputsuffix),sprintf('sq-%d%s.set',subject,outputsuffix),duds);

clear *suffix duds

%% EPOCH DATA

Constants = GetSquaresConstants;
handles.subject = subject;
handles.version = 'noduds';
handles.eventnumbers = Constants.SACCADEEND_BASE+[Constants.DISTRACTOR, Constants.INTEGRATION, Constants.COMPLETION];
handles.eventnames = {'distractor','integration','completion'};
handles.reference = 'None';
handles.epochtimes = [-1000 2000];
handles.baselinetimes = [-200 0];
handles.baseevent = 'Each Event';
handles.removeIncorrect = false;
handles.removeBlinks = false;
handles.useCutoff = false;
handles.tepoch = [-500 500];

filename = sprintf('sq-%d-all-%s.set',handles.subject,handles.version); % use suffix indicated by input string 'version'

% USE BASELINE RELATIVE TO STIMULUS APPEARANCE
if strcmp(handles.baseevent,'Stimulus')
    % epoch to appear time with baseline subtraction and wide epoch times
    [ALLEEG] = EpochEeglabFile(filename,[-2000 4000],handles.baselinetimes,...
        Constants.TRIAL_START,{'baseline=appear'},handles.reference); 
    % re-epoch to individual events with desired epoch times and no baseline subtraction
    [ALLEEG, EEG, CURRENTSET] = EpochEeglabFile(ALLEEG(end),handles.epochtimes,[],...
        handles.eventnumbers,handles.eventnames,handles.reference); 

% OR USE BASELINE RELATIVE TO EACH INDIVIDUAL EVENT
else
    [ALLEEG, EEG, CURRENTSET] = EpochEeglabFile(filename,handles.epochtimes,handles.baselinetimes,handles.eventnumbers,handles.eventnames,handles.reference);
end
    
% PERFORM EPOCHING
nFiles = length(ALLEEG);
for i=2:nFiles
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',i,'study',0);
    if handles.removeIncorrect
        EEG = RemoveIncorrectTrials(EEG);
    end
    if handles.removeBlinks
        EEG = RemoveEyeBlinkTrials(EEG,Constants.BLINKSTART,handles.tepoch,handles.tsaccade,handles.tdiscrim);        
    end
    if handles.useCutoff
        EEG = ExcludeTrialsByCutoff(EEG,handles.cutoff,handles.tepoch,handles.tsaccade,handles.tdiscrim);
    end
    %  Store results (but do not save)
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off'); 
    EEG = eeg_checkset( EEG );
end
% Send results to base workspace so user can explore (and use them for
% analysis parts of this gui)
% assignin('base','ALLEEG',ALLEEG);
% assignin('base','EEG',EEG);
% assignin('base','CURRENTSET',CURRENTSET);
% evalin('base','eeglab redraw');
eeglab redraw