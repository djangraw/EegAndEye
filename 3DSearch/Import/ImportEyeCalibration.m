function [EEG, output_filename]  = ImportEyeCalibration(subject,session,filetype)

% Imports an EEG data file produced by the eyecalibration.es script
% (from the LIINC website) into EEGLAB, adds electrode locations, and adds
% the events in a way that works for the pop_eyesubtract plugin later.
% 
% [EEG, output_filename] = ImportEyeCalibration(subject,session,filetype)
%
% - EEGLAB should already be started.
% 
% INPUTS:
% -subject and session are scalars referring to the eyecalibration file
% that was recorded.
% -filetype is a string indicating the eeg system used to record that file.
%   Sensorium: eyecalibration-<subject>-<session>.b.dat
%   Biosemi: eyecalibration-<subject>-<session>.bdf
%
% OUTPUTS:
% -EEG is the eeglab struct containing the imported eyecalibration data.
% -output_filename is a string indicating where that data was saved (for
% easy input to UseEyeCalibration() ).
%
% Created 8/20/10 by DJ.
% Updated 11/23/10 by DJ - added filetype switch for Sensorium data
% Updated 11/30/10 by DJ - fixed Sensorium scaling issue, removed
% re-referencing to average
% Updated 12/2/10 by DJ - added filtering and downsampling!
% Updated 12/16/10 by DJ - cleanup
% Updated 2/21/11 by DJ - adjusted to new rereferenceData input. 
% Updated 2/28/11 by DJ - rereferenceData is now a function
% Updated 3/2/11 by DJ - made function, moved reref to UseEyeCalibration().
% Updated 7/28/11 by DJ - added option for new Sensorium cap (filetype =
% 'Sensorium-2011').

%% CHECK INPUTS AND SET UP
% Handle defaults
if nargin<3
    filetype = 'Biosemi';
end


data_dir = [cd '/'];
% data_dir = '/Users/dave/Documents/Data/3DSearch/';
electrodeloc_dir = '/Users/dave/Documents/Tools/eeglab/locationfiles/';

% Set filter bounds
highpass_cutoff = 0.5;
notch_bounds = [59 61];
lowpass_cutoff = 100;

% refChans = 'None'; % use this line to skip re-referencing
% refChans = ''; % use this line to re-reference to average

% Clear datasets from EEGLAB
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
eeglab redraw;

fprintf('---Importing eyecalibration data for subject %d, session %d---\n',subject,session);
tic


%% Load the data
switch filetype
    case 'Biosemi'
        % Set up
        electrodeloc_file = 'biosemi64.sph';
        % Import data
        EEG = pop_biosig(sprintf('%seyecalibration-%d-%d.bdf',data_dir,subject,session));
        EEG.setname = sprintf('eyecalibration-%d-%d',subject,session);
        EEG = eeg_checkset( EEG );[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',sprintf('eyecalibration-%d-%d',subject,session),'gui','off'); 
        % pop_eegplot( EEG, 1, 1, 1); % view the raw data
        
        % Get port events
        dat = sopen(sprintf('%seyecalibration-%d-%d.bdf',data_dir,subject,session));
        port_events = [dat.BDF.Trigger.POS, dat.BDF.Trigger.TYP-16128, [diff(dat.BDF.Trigger.POS); 0]];
        raw_fs = dat.SampleRate;
        new_fs = 256;
        
        % clean up
        clear dat;
    case 'Sensorium-2008'
        
        %Set up
        electrodeloc_file = 'jen_sensorium87_79chan.loc';
        input_filename = sprintf('eyecalibration-%d.%.3d.b',subject,session);
        channels = [1:60 67:70 72 74:87]; % exclude EOGs, A1, A2 (see 'Sajda87 ch capfinal.xls')
        raw_fs = 1000;
        new_fs = 250;
        nSamples = []; % autodetect length of file
        offset = 0;
        scale = 10^6;
        
        % Read in data
        [eegdata, events] = readEEG_b4preprocessing(input_filename,channels,nSamples,offset);
        % Scale to microvolts
        eegdata = eegdata*scale; % convert from V to uV
        % Import
        EEG = pop_importdata('dataformat','array','nbchan',0,'setname',sprintf('eyecalibration-%d-%d',subject,session),...
            'srate',raw_fs,'subject',num2str(subject),'pnts',0,'xmin',0,'session',num2str(session));
        EEG.data = eegdata;
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
        EEG = eeg_checkset( EEG );
        
        % Get port events
        % arrange events into a matrix where the columns are [latency, type]
        iEvents = find(diff(events)~=0)+1;
        durations = [diff(iEvents) 0];
        port_events = [iEvents', events(iEvents)' durations'];
        
        % clean up
        clear input_filename channels nSamples offset scale iEvents times durations eegdata events
    case 'Sensorium-2011'
        
        %Set up
        electrodeloc_file = 'dave_sensorium90_84chan.loc';
        input_filename = sprintf('eyecalibration-%d.%.3d.b',subject,session);
        channels = [7:90]; % exclude EOGs (see 'Sajda10-10 cap.xlsx')
        raw_fs = 1000;
        new_fs = 250;
        nSamples = []; % autodetect length of file
        offset = 0;
        scale = 10^6;
        
        % Read in data
        [eegdata, events] = readEEG_b4preprocessing(input_filename,channels,nSamples,offset);
        % Scale to microvolts
        eegdata = eegdata*scale; % convert from V to uV
        % Import
        EEG = pop_importdata('dataformat','array','nbchan',0,'setname',sprintf('eyecalibration-%d-%d',subject,session),...
            'srate',raw_fs,'subject',num2str(subject),'pnts',0,'xmin',0,'session',num2str(session));
        EEG.data = eegdata;
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
        EEG = eeg_checkset( EEG );
        
        % Get port events
        % arrange events into a matrix where the columns are [latency, type]
        iEvents = find(diff(events)~=0)+1;
        durations = [diff(iEvents) 0];
        port_events = [iEvents', events(iEvents)' durations'];
        
        % clean up
        clear input_filename channels nSamples offset scale iEvents times durations eegdata events
    otherwise
        error('filetype %s not recognized!',filetype);
end





%% Put in the channel location information
EEG=pop_chanedit(EEG, 'load',{[electrodeloc_dir, electrodeloc_file] 'filetype' 'autodetect'});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

%% Low-Pass Filter the data
[EEG, ~, lpf] = pop_eegfilt_returnfilter( EEG, 0, lowpass_cutoff, [], [0]); %low-pass, 100Hz cutoff
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',sprintf('3DS-%d-%d-lowpass',subject,session),'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );

%% Resample data to new_fs Hz
EEG = pop_resample_myfilter( EEG, new_fs, lpf);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',sprintf('3DS-%d-%d-resampled',subject,session),'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );

%% Filter the data further
EEG = pop_eegfilt( EEG, highpass_cutoff, 0, [], [0]); %high-pass, 0.5Hz cutoff
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',sprintf('3DS-%d-%d-highpass',subject,session),'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_eegfilt( EEG, notch_bounds(1), notch_bounds(2), [], [1], [1]); %notch, 59-61Hz, FFT
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',sprintf('3DS-%d-%d-filtered',subject,session),'overwrite','on','gui','off'); 

%% Add events
assignin('base','port_events',port_events);
EEG = pop_importevent( EEG, 'append','no','event','port_events','fields',{'latency' 'type' 'duration'},'timeunit',1/raw_fs,'optimalign','off');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
evalin('base','clear port_events');

%% Save results
output_filename = sprintf('eyecalibration-%d-%d.set',subject,session);
EEG = pop_saveset( EEG, 'filename',output_filename,'filepath',data_dir);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% Clean up
toc % Display elapsed time
eeglab redraw
clear eegdata events data_dir electrodeloc_* dat port_events raw_fs new_fs highpass_cutoff lowpass_cutoff notch_bounds lpf

