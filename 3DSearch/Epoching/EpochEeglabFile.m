function [ALLEEG, EEG, CURRENTSET] = EpochEeglabFile(filename,epoch_times,baseline_times,event_numbers,event_names,refChans)

% Takes an eeglab file with events attached, and creates a new dataset
% consisting of epochs around each event type specified.
%
% [ALLEEG, EEG, CURRENTSET] = EpochEeglabFile(filename,epoch_times,
%     baseline_times,event_numbers,event_names,refChans)
%
% - An example event is Numbers.SACCADE_TO+Numbers.TARGET (41).
% - EEGLAB should already be started so we have access to their functions.
%
% INPUTS:
% - filename is a string specifying a .set file in the current path, or an
% EEGLAB struct that has already been loaded.  If an array of structs is
% given, the last one will be used for epoching, but all will remain in the
% ALLEEG output struct.
% - epoch_times is a 2-element vector specifying the start and end times 
% (ms) that each epoch should extend around the anchoring event [-500 1000]
% - baseline_times is a 2-element vector specifying the start and end times
% (ms) that each epoch should use for baseline subtraction.  If
% baseline_times is empty, no baseline subtraction will be performed.
% - event_numbers is an n-element vector of event codes that will be used
% to anchor the epochs.  They should match the codes in the eeglab file.
% - event_names is an n-element cell array of strings, each of which gives
% a short text tag for the corresponding event_number. If event_names is a
% single cell containing a string, all events will be used at the same time
% to create a single epoched output struct.
% - refChans may be used to define which channels will be used
% to rereference the data (see rereferenceData.m for details).
%
% OUTPUTS:
% -ALLEEG, EEG, and CURRENTSET are the standard eeglab variables of the
% same names.
%
% Created 8/3/10 by DJ.
% Updated 8/19/10 by DJ - added Numbers constants, removed button presses.
% Updated 8/20/10 by DJ - comments.
% Updated 11/2/10 by DJ - added re-referencing to average
% Updated 11/4/10 by DJ - rereferenceData, not rereferenceToAvg
% Updated 11/8/10 by DJ - using old_setname
% Updated 11/23/10 by DJ - switched baseline to [-200 0]
% Updated 12/17/10 by DJ - switched epoch_times to [-0.5 1]
% Updated 2/21/11 by DJ - adjusted to new rereferenceData input, added
%  'version' input. 
% Updated 2/28/11 by DJ - made function. (rereferenceData also a fn now).
% Updated 3/1/11 by DJ - comments.
% Updated 3/28/11 by DJ - first input can be an eeglab struct, empty 
%  baseline_times means no baseline subtraction, single-cell event_names
%  means single-struct output with all events treated equally.


%% CHECK INPUTS AND SET UP

GetNumbers;
if nargin<2 || isempty(epoch_times)
    epoch_times = [-500 1000]; % time, in ms, that each epoch should extend around the anchoring event
end
if nargin<3
    baseline_times = [-200 0]; % time, in ms, relative to the anchoring event, that should be considered baseline.
end
if nargin<4 
    event_numbers = [Numbers.SACCADE_TO+Numbers.TARGET, Numbers.SACCADE_TO+Numbers.DISTRACTOR,...
        Numbers.ENTERS+Numbers.TARGET, Numbers.ENTERS+Numbers.DISTRACTOR,...
        Numbers.SEES+Numbers.TARGET, Numbers.SEES+Numbers.DISTRACTOR,...
        Numbers.BUTTON, Numbers.SIM_BUTTON]; %See script GetNumbers for values of these codes
end
if nargin<5 
    event_names = {'targsac','distsac','targapp','distapp','targsee',...
        'distsee','targbutton','disbutton'}; % used to name new datasets
end
if nargin<6 
    refChans = 'None';
%     refChans = ''; % use this line to re-reference to average
end

if numel(event_names)==1
    do_singleoutput = true;
else
    do_singleoutput = false;
end

% Data info
data_dir = [cd '/'];
% data_dir = '/Users/dave/Documents/Data/3DSearch/';

% Clear datasets from EEGLAB
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
eeglab redraw;

tic

if isnumeric(event_numbers)
    for i=1:numel(event_numbers)
        event_strings{i} = num2str(event_numbers(i));
    end
else
    event_strings = event_numbers;
end
%% Load dataset
if ischar(filename)
    fprintf('---Epoching file %s into %d datasets around events %s ---\n',filename,numel(event_names),sprintf('%s ',event_strings{:}));
    EEG = pop_loadset('filename',filename,'filepath',data_dir);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 ); % add to ALLEEG
elseif isstruct(filename)
    ALLEEG = filename;
    EEG = ALLEEG(end);
    fprintf('---Epoching file %s into %d datasets around events %s ---\n',EEG.setname,numel(event_names),sprintf('%s ',event_strings{:}));
else
    error('first input must be either string or struct!');
end

%% Re-reference to specified channels
EEG = rereferenceData(EEG,refChans);
old_setname = EEG.setname;

%% Epoch data
if do_singleoutput
    % create a single new dataset in which all the epoch events are treated equally
    % convert event_numbers to a cell array of strings
    EEG = pop_epoch( EEG, event_strings, epoch_times/1000, 'newname', sprintf('%s-%s',old_setname,event_names{1}), 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
    EEG = eeg_checkset( EEG );
    if ~isempty(baseline_times) % if baseline_times is empty, skip baseline subtraction!
        EEG = pop_rmbase( EEG, baseline_times); % otherwise, subtract the baseline times indicated by user.
    end
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, numel(event_strings)+1,'retrieve',1,'study',0); 
    EEG = eeg_checkset( EEG );
else
    % Create a new dataset for each epoch group
    for i=1:numel(event_numbers)
        EEG = pop_epoch( EEG, event_strings(i), epoch_times/1000, 'newname', sprintf('%s-%s',old_setname,event_names{i}), 'epochinfo', 'yes');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
        EEG = eeg_checkset( EEG );
        if ~isempty(baseline_times) % if baseline_times is empty, skip baseline subtraction!
            EEG = pop_rmbase( EEG, baseline_times); % otherwise, subtract the baseline times indicated by user.
        end
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, i+1,'retrieve',1,'study',0); 
        EEG = eeg_checkset( EEG );
    end
end
%% Clean up
toc % Display elapsed time

