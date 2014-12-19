% RunFlightSimPilot_script.m
%
% Runs through the steps of importing and analyzing Flight Simulator/PIO
% project pilot data.
%
% Created 9/23/14 by DJ.
% Updated 10/22/14 by DJ - added reject-session section.

% Declare parameters for this experiment/subject
experiment = 'FlightSim_1pt0';
subject = 3;
switch subject
    case 1
        sessions = 1:50;
        eegSessions = [1:7 16 26 36 46];
    case {2,3}
        sessions = 1:40; % NEDE sessions
        eegSessions = [1 11 21 31]; % BDF filenames
    otherwise
        disp('Please enter this subject''s sessions into the script.')
        sessions = [];
        eegSessions = [];
end

% ------------------------ %
% ----- BEHAVIOR DATA ---- %
% ------------------------ %

%% Import behavior
% Save behavior .mat files
NEDE_ImportData(experiment,subject,sessions);

%% Import Raw Pupil/Gaze Samples
% Save eye position/pupil size samples
cd samples
GetEyeSamples(subject,sessions,experiment,0);
cd ..

%% Load behavior results
% Load into a vector of NEDE structs
clear y
for i=1:numel(sessions);
    load(sprintf('%s-%d-%d',experiment,subject,sessions(i)));
    y(i) = x;
end
save(sprintf('%s-%d-all',experiment,subject),'y');

% Make combo struct with all trials appended but only 1 copy of params/objects fields
eventsstruct = AppendStructs({y.events},1); % only events field will be different across files
combostruct = y(1); % keep only 1 copy of other fields
combostruct.events = eventsstruct;

%% Get pupil data and create flight sim behavior file
nTrials = length(x.events.trial.time_start);
pup_cell = cell(1,nTrials);
for i=1:numel(sessions);    
    load(sprintf('%s-%d-%d-eyepos.mat',experiment,subject,sessions(i)));
    pup_cell{i} = pupilsize;
end

filename = GetFlightSimBehavior(combostruct,pup_cell); % save results as <expt>-<subj>-all-Behavior.mat

%% Plot Behavior results

% filename = sprintf('%s-%d-all-Behavior.mat',experiment,subject);
PlotFlightSimBehavior(filename); % trajectories, PIO amplitude, pupil size plotted for all trials superimposed

%% Get/Plot pre-PIO epochs
pioThreshold = 100; % (arbitrary) threshold for Hilbert transform amplitude at which subject enters a PIO
PlotPrePioBehavior(filename,pioThreshold); % grayscale histograms with red mean line on top


% ------------------------ %
% ------- EEG DATA ------- %
% ------------------------ %
%% Import EEG data
% set filter and sampling parameters
bandpass_bounds = [1 100];
notch_bounds = [59 61];
new_fs = 256;

% Import and filter data.
for i=1:numel(eegSessions)
    % Import and save in eeglab
    NEDE_ImportToEeglab(experiment,subject,eegSessions(i),'Biosemi','-filtered',bandpass_bounds,notch_bounds,new_fs);
end
% Combine EEG data and save result
CombineEeglabSessions(subject,eegSessions,'-filtered','-all-filtered',experiment);

%% Reject bad sessions
switch subject
    case 2
        badSessions = [18 21 37];
    case 3
        badSessions = [2 7 9 11];
    otherwise
        badSessions = [];
end
% check whether bad sessions have already been deleted
load(sprintf('%s-%d-all',experiment,subject));
if numel(y)<numel(sessions)
    disp('Some sessions have already been deleted! Not deleting any more.')
else
    % remove from behavior struct array
    fprintf('Deleting %d sessions from behavior...\n',numel(badSessions));
    y(badSessions) = [];
    save(sprintf('%s-%d-all',experiment,subject),'y');
    % remove from EEG struct
    fprintf('Deleting %d sessions from EEG...\n',numel(badSessions));
    EEG = pop_loadset(sprintf('%s-%d-all-filtered.set',experiment,subject));
    EEG = NEDE_RejectSession(EEG, badSessions);
    EEG = pop_saveset(EEG,'filename',sprintf('%s-%d-all-filtered.set',experiment,subject));
end
disp('Done!')

%% Load combined file
% load the file we just saved
EEG = pop_loadset(sprintf('%s-%d-all-filtered.set',experiment,subject));

%% Sync up files
% Declare params
newEegFilename = sprintf('%s-%d-all-events.set',experiment,subject);
SYNC_CODES = [211, 208]; % include 208 for subj 1 (last 2 bits of event codes lost in EEG file)
% Check how well codes synched up
NEDE_CheckSync(y,EEG,SYNC_CODES);
% Import standard NEDE events like saccades, blinks, etc.
[eyeTimes, eyeCodes] = NEDE_GetEeglabEvents(y,EEG,SYNC_CODES);
EEG = NEDE_AddEeglabEvents(y,EEG,eyeTimes,eyeCodes,SYNC_CODES);
% Save dataset with events
pop_saveset(EEG,'filename',newEegFilename); 


%% Add PIO Events
% THIS CODE SHOULD BE REPLACED WHEN A MORE NUANCED DEFINITION OF PIO IS
% IMPLEMENTED!!!!

% Declare threshold on Hilbert Transform amplitude for entering PIO
pioThreshold = 100;
% Load ring positions from behavior file
foo = load(sprintf('%s-%d-all-Behavior',experiment,subject),'ringpos');
ringpos = foo.ringpos;

% Calculate times of PIO onsets
[eyeTimes,eyeCodes] = deal(cell(1,numel(y)));
for i=1:numel(y)
    % Get camera positions & times
    zpos = y(i).events.camera.position(:,2);
    ypos = y(i).events.camera.elevation;
    tim = y(i).events.camera.time;
    % Get ideal path from ring position record
    idealy = interp1(ringpos(:,1),ringpos(:,2),zpos,'linear','extrap');    
    % Get amplitude of PIO
    patherror = abs(ypos-idealy);
    patherror(zpos<0) = 0;
    hilbpatherror = abs(hilbert(patherror));
    hilbpatherror(zpos<0) = 0; % Avoid edge effects in hilbert transform (takeoff had no error)
    % Find first time when pio amplitude exceeded threshold
    iPio = find(hilbpatherror>pioThreshold,1);
    % if we never exceeded the threshold, use the last sample as the PIO
    % onset time.
    if isempty(iPio)
        iPio = length(hilbpatherror);
    end
    % Log in times/codes cells for input to NEDE_AddEeglabEvents
    eyeTimes{i} = tim(iPio);
    eyeCodes{i} = repmat({'StartPio'},size(eyeTimes{i}));
end

SYNC_CODES = [211, 208]; % include 208 for subj 1 (last 2 bits of event codes lost in EEG file)
EEG = pop_loadset(sprintf('%s-%d-all-events.set',experiment,subject)); % load eeglab set
EEG = NEDE_AddEeglabEvents(y,EEG,eyeTimes,eyeCodes,SYNC_CODES); % add events (no saving)



%% Epoch around PIO events
% Declare epoch parameters
tEpoch_ms = [-5000 0]; % start and end time of epoch, in ms
tBaseline_ms = [-100 0]; %start and end time of baseline window, in ms
% Epoch data
EEG_epoch = pop_epoch( EEG, {  'StartPio'  }, tEpoch_ms/1000, 'newname', 'FlightSim_1pt0-1-all epochs', 'epochinfo', 'yes');
% Remove baseline
EEG_epoch = pop_rmbase( EEG_epoch, tBaseline_ms);

%% Plot scalp maps
% Declare plot parameters
tPlots_ms = -475:50:-25; % times of bin centers of plots (in ms)
tWidth_ms = 50; % width of time bins (in ms)

% extract averages within time bins
[~,sm_all] = GetScalpMaps(EEG_epoch.data,EEG_epoch.times,tPlots_ms,tWidth_ms);
% Plot scalp maps
PlotScalpMaps(sm_all,EEG_epoch.chanlocs,[],tPlots_ms);