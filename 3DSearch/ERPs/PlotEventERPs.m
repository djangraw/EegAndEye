% PlotEventERPs
%
% Created 10/28/10 by DJ (as PlotBlinkERPs.m).
% Updated 11/2/10 by DJ - renamed PlotEventERPs.m, added event type switch

% eventType = 'blink';
% eventType = 'target';
eventType = 'distractor';

%% CHECK INPUTS AND SET UP
if ~exist('subject','var')
    error('variable ''subject'' must be defined to use this script!');
elseif ~exist('ALLEEG','var')
    error('Start EEGLAB before using this script!');
end

% Data info
data_dir = [cd '/'];
% data_dir = '/Users/dave/Documents/Data/3DSearch/';
filename = sprintf('3DS-%d-all-filtered.set',subject); % combined dataset as made in CombineEeglabSessions
% filename = sprintf('3DS-%d-all-eyeremoved.set',subject); % eyecalibration-removed dataset as made in UseEyeCalibration
% filename = sprintf('3DS-%d-all-nofrontal.set',subject); % removed frontal electrodes

% Event info
GetNumbers;
epoch_times = [-1 2]; % time, in seconds, that each epoch should extend around the anchoring event
baseline_times = [-300 0]; % time, in ms, relative to the anchoring event, that should be considered baseline.
switch eventType
    case 'blink'
        event_number = [Numbers.BLINK]; %See script GetNumbers for values of these codes
        event_name = 'blink'; % used to name new datasets

        % ERP info
        plot_times = [-300 500]; % time range, in ms, relative to the anchoring event, that will get an ERP.
        scalpmap_times = [0 100 200]; % time, in ms, relative to the anchoring event, that will get scalp maps.
        channel_toplot = 'Cz'; % the channel number or name that gets single-trial plotting
    case 'target'
        event_number = [Numbers.TARGET+Numbers.ENTERS]; %See script GetNumbers for values of these codes
        event_name = 'targetAppear'; % used to name new datasets

        % ERP info
        plot_times = [-500 1500]; % time range, in ms, relative to the anchoring event, that will get an ERP.
        scalpmap_times = [300:100:600]; % time, in ms, relative to the anchoring event, that will get scalp maps.
        channel_toplot = 'Pz'; % the channel number or name that gets single-trial plotting
    case 'distractor'
        event_number = [Numbers.DISTRACTOR+Numbers.ENTERS]; %See script GetNumbers for values of these codes
        event_name = 'distAppear'; % used to name new datasets

        % ERP info
        plot_times = [-500 1500]; % time range, in ms, relative to the anchoring event, that will get an ERP.
        scalpmap_times = [300:100:600]; % time, in ms, relative to the anchoring event, that will get scalp maps.
        channel_toplot = 'Pz'; % the channel number or name that gets single-trial plotting
end

% Clear datasets from EEGLAB
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
eeglab redraw;

fprintf('---Getting blink info from file %s ---\n');
tic

%% Load dataset
EEG = pop_loadset('filename',filename,'filepath',data_dir);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% Get channel number
if isstr(channel_toplot)
    channel_toplot = find(strcmp(channel_toplot,{EEG.chanlocs.labels}));
end

%% Re-reference to average
rereferenceToAvg;

%% Epoch data
% Create a new dataset
EEG = pop_epoch( EEG, { num2str(event_number) }, epoch_times, 'newname', sprintf('3DS-%d-all-%s',subject,event_name), 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, baseline_times);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% Plot blink ERP for all channels
figure;
pop_timtopo(EEG, plot_times, scalpmap_times, sprintf('ERP data and scalp maps of %s',filename));

%% Plot blink ERP and single-trial sutff for one channel
figure; 
pop_erpimage(EEG,1, channel_toplot,[],EEG.chanlocs(channel_toplot).labels,1,1,{},[],'' ,'yerplabel','\muV','erp','on','limits',[plot_times NaN NaN NaN NaN NaN NaN] ,'cbar','on','topo', { channel_toplot EEG.chanlocs EEG.chaninfo } );

%% Clean up
toc % Display elapsed time
eeglab redraw;
clear data_dir filename epoch_times baseline_times plot_times event_number event_name Numbers



