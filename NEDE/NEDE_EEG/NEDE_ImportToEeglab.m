function NEDE_ImportToEeglab(prefix,subject,eegSession,eegFiletype,output_suffix,bandpass_bounds,notch_bounds,new_fs,include_eog)

% NEDE_ImportToEeglab(prefix,subject,eegSession,eegFiletype,output_suffix,
%  bandpass_bounds,notch_bounds,new_fs,include_eog)
% 
%
% -Imports an eeg data file from a 3DSearch experiment into EEGLAB, then
% filters and downsamples it to make a more compact file that's ready for 
% analysis.
% -EEGLAB should already be started when you run this script.
%
% INPUTS:
% - subject and eegSession are numbers identifying the data file. They
% must be defined in the workspace when this script is called.
% - eegFiletype is a string indicating the type of eeg file being used as
% input. {'Biosemi' 'Sensorium-2008' 'Sensorium-2011'}
% - prefix is a string indicating which experiment is being used. This will
% influence the input and output filenames.
% - output_suffix is a string indicating the end of the output data file
% (see outputs section below). Default: '-filtered'
% - bandpass_bounds is a 2-element vector indicating the high-pass and
% low-pass cutoffs (in Hz) of the band-pass filter. Default = [0.5, 100]
% - notch_bounds is a 2-element vector indicating the low and high cutoffs
% (in Hz) of the notch filter. Default = [59, 61]
% - new_fs is a scalar indicating the sampling rate you'd like to
% downsample to.  Default = 256 for Biosemi, 250 for Sensorium.
%
% OUTPUTS:
% - The result will be saved in the directory data_dir (default is 
% current directory). Edit code directly to change.
%
% Created 7/29/10 by DJ as ImportToEeglab.
% Updated 8/20/10 by DJ - comments
% Updated 8/31/10 by DJ - resolved problem of unwanted filtering during
%   pop_resample by implementing custom functions pop_eegfilt_returnfilter
%   and pop_resample_myfilter.
% Updated 10/15/10 by DJ - moved re-referencing to after filtering
% Updated 11/2/10 by DJ - removed re-referencing!
% Updated 11/30/10 by DJ - made into a function, added Sensorium support
% Updated 12/15/10 by DJ - added ascSession/eegSession
% Updated 7/28/11 by DJ - added new Sensorium cap option (eegFileType =
%  'Sensorium-2011')
% Updated 10/27/11 by DJ - added experimentType input to work with Squares
%  experiments
% Updated 12/19/11 by DJ - added filter bounds inputs
% Updated 1/6/12 by DJ - added Sensorium-2011-50Hz option
% Updated 6/29/12 by DJ - added Sensorium-2011-1000Hz option & allow no lpf/hpf
% Updated 7/3/12 by DJ - added new_fs and include_eog as inputs, set up iir
%  filter iff. low-pass cutoff <= 100Hz
% Updated 2/28/13 by DJ - added SquaresFix compatibility, Sensorium-2013.
% Updated 5/17/13 by DJ - added SquaresFix3 compatibility.
% Updated 2/25/14 by DJ - switched to sensorium loc files that include PO5.
% Updated 9/25/14 by DJ - switched name to NEDE_...
% Updated 10/17/14 by DJ - fixed bug where biosemi events are offset by 16128.
% Updated 10/22/14 by DJ - used NEDE_ExtractBiosemiEvents

%% CHECK INPUTS AND SET UP

if ~exist('output_suffix','var') || isempty(output_suffix)
    output_suffix = '-filtered';
end
if ~exist('bandpass_bounds','var')  || isempty(bandpass_bounds)
    bandpass_bounds = [0.5 100];
end
if ~exist('notch_bounds','var')  || isempty(notch_bounds)
    notch_bounds = [59 61];
end
if ~exist('new_fs','var') 
    new_fs = [];
end
if ~exist('include_eog','var') 
    include_eog = 0;
end

% Set filter bounds
highpass_cutoff = bandpass_bounds(1);
lowpass_cutoff = bandpass_bounds(2);

% make sure eeglab is running
eeglab_is_running = evalin('base','exist(''ALLEEG'',''var'')'); % if ALLEEG is a variable in the base class, eeglab is running.
if ~eeglab_is_running
    error('Start EEGLAB before using this function!');
end

% Directories
data_dir = [cd '/'];
% data_dir = '/Users/dave/Documents/Data/3DSearch/';
% electrodeloc_dir = '/Users/dave/Documents/Tools/eeglab/locationfiles/'; % on Tesla machine
electrodeloc_dir = '/Users/jangrawdc/Documents/MATLAB/LIINC/locationfiles/'; % on DJ's NIH machine

% Filenames
filename_base = sprintf('%s-%d-%d',prefix,subject,eegSession);
final_filename = sprintf('%s%s.set',filename_base,output_suffix); 


% Set EEG-system-specific constants
switch eegFiletype
    case 'Biosemi'
        eeg_filename = sprintf('%s%s-%d-%d.bdf',data_dir,prefix,subject,eegSession);
        electrodeloc_filename = 'biosemi64.sph';
        channels = 1:64;        
        raw_fs = 2048; % sampling frequency of raw data
        if isempty(new_fs)
            new_fs = 256;
        end
        lpf_useIIR = 0;
        lpf_order = [];
    case 'Sensorium-2008'
        eeg_filename = sprintf('%s-%d.%.3d.b',prefix,subject,eegSession);
        electrodeloc_filename = 'jen_sensorium87_79chan.loc';
        channels = [1:60 67:70 72 74:87]; % exclude EOGs, A1, A2 (see 'Sajda87 ch capfinal.xls')
        nSamples = []; % autodetect length of file
        offset = 0;
        scale = 10^6; % output of readEEG_b4preprocessing is in V, eeglab expects uV
        raw_fs = 1000; % sampling frequency of raw data
        if isempty(new_fs)
            new_fs = 250;
        end
        lpf_useIIR = 0;
        lpf_order = [];
    case {'Sensorium-2011' 'Sensorium-2013'}
        if strcmp(eegFiletype,'Sensorium-2011')
            eeg_filename = sprintf('%s-%d.%.3d.b',prefix,subject,eegSession);
            if include_eog
%                 channels = 1:90;
                channels = 1:91;
            else
%                 channels = 7:90;
                channels = 7:91;
            end
        else
            eeg_filename = sprintf('%s-%d.%.3d.bh',prefix,subject,eegSession);
            if include_eog
%                 channels = 3:92;
                channels = 3:93;
            else
%                 channels = 9:92;
                channels = 9:93;
            end
        end    
        if include_eog
%             electrodeloc_filename = 'dave_sensorium90chan.loc';
            electrodeloc_filename = 'dave_sensorium91chan.loc';
        else
%             electrodeloc_filename = 'dave_sensorium90_84chan.loc';
            electrodeloc_filename = 'dave_sensorium91_85chan.loc';
        end
        nSamples = []; % autodetect length of file
        offset = 0;
        scale = 10^6; % output of readEEG_b4preprocessing is in V, eeglab expects uV
        raw_fs = 1000; % sampling frequency of raw data        
        if isempty(new_fs)
            new_fs = 250;   % 1000/6;  % 1000; 
        end
        if bandpass_bounds(2) <= 100 % low low-pass cutoff --> sharp filter
            disp('Using 12th-order iir filter...')
            lpf_useIIR = 1;
            lpf_order = 12;
        else
            disp('Using fir filter of unspecified order...')
            lpf_useIIR = 0;
            lpf_order = [];
        end
        
    % SPECIAL PRESETS
    case 'Sensorium-2011-50Hz'
        eeg_filename = sprintf('%s-%d.%.3d.b',prefix,subject,eegSession);
%         electrodeloc_filename = 'dave_sensorium90_84chan.loc';
%         channels = [7:90]; % exclude EOGs (see 'Sajda10-10 cap.xlsx')
        electrodeloc_filename = 'dave_sensorium91_85chan.loc';
        channels = 7:91; % exclude EOGs (see 'Sajda10-10 cap.xlsx')
        nSamples = []; % autodetect length of file
        offset = 0;
        scale = 10^6; % output of readEEG_b4preprocessing is in V, eeglab expects uV
        raw_fs = 1000; % sampling frequency of raw data
%         new_fs = raw_fs/floor(raw_fs/(2*lowpass_cutoff)); 
        if isempty(new_fs)
            new_fs = 1000/6;
        end
        lpf_useIIR = 1;
        lpf_order = 12;
    case 'Sensorium-2011-1000Hz'
        eeg_filename = sprintf('%s-%d.%.3d.b',prefix,subject,eegSession);
%         electrodeloc_filename = 'dave_sensorium90_84chan.loc';
%         channels = [7:90]; % exclude EOGs (see 'Sajda10-10 cap.xlsx')
        electrodeloc_filename = 'dave_sensorium91_85chan.loc';
        channels = 7:91; % exclude EOGs (see 'Sajda10-10 cap.xlsx')
        nSamples = []; % autodetect length of file
        offset = 0;
        scale = 10^6; % output of readEEG_b4preprocessing is in V, eeglab expects uV
        raw_fs = 1000; % sampling frequency of raw data   
        if isempty(new_fs)
            new_fs = 1000;      
        end
        lpf_useIIR = 0;
        lpf_order = [];
        
    otherwise
        error('FileType not recognized!');
end

tic;
fprintf('---Importing data for subject %d, session %d and saving to %s---\n',subject,eegSession,final_filename);
ALLEEG = [];

%% Load the data
switch eegFiletype
    case 'Biosemi'
        % Import (EEGLAB does the work for us)
        EEG = pop_biosig(eeg_filename,'channels',channels,'ref',channels);
        % write events to temporary file for reading into EEGLAB
        events_filename = 'TEMP_events.txt';
        events = NEDE_ExtractBiosemiEvents(eeg_filename,events_filename);
        % import events to eeglab
        EEG = pop_importevent( EEG, 'append','no','event',events_filename,'fields',{'latency' 'type'},'timeunit',1/EEG.srate,'optimalign','off');
        delete(events_filename); % remove temporary file

    case {'Sensorium-2008' 'Sensorium-2011' 'Sensorium-2011-1000Hz' 'Sensorium-2011-50Hz'}
        % Get data
        eegdata = readEEG_b4preprocessing(eeg_filename,channels,nSamples,offset); % An's program for importing Sensorium files
        % Scale to microvolts
        eegdata = eegdata*scale; % convert from V to uV
        % Import
        EEG = pop_importdata('dataformat','array','nbchan',0,'srate',raw_fs,...
            'subject',num2str(subject),'pnts',0,'xmin',0,'session',num2str(ascSession));
        EEG.data = eegdata;
    case 'Sensorium-2013'
        % Get data
        eegdata = readEEG_sensorium2013(eeg_filename,channels); % simple program for importing Sensorium files
        % Scale to microvolts
        eegdata = eegdata*scale; % convert from V to uV
        % Import
        EEG = pop_importdata('dataformat','array','nbchan',0,'srate',raw_fs,...
            'subject',num2str(subject),'pnts',0,'xmin',0,'session',num2str(ascSession));
        EEG.data = eegdata;
    otherwise
        error('FileType not recognized!');
end
EEG.setname = sprintf('%s-%d-%d-raw',prefix,subject,eegSession);
EEG = eeg_checkset( EEG );
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',[filename_base '-raw'],'gui','off'); 
% pop_eegplot( EEG, 1, 1, 1); % view the raw data


%% Put in the channel location information
EEG=pop_chanedit(EEG, 'load',{[electrodeloc_dir electrodeloc_filename] 'filetype' 'autodetect'});
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

%% Low-Pass Filter the data
if ~isnan(lowpass_cutoff) && highpass_cutoff < lowpass_cutoff
    [EEG, ~, lpf] = pop_eegfilt_returnfilter( EEG, 0, lowpass_cutoff, lpf_order, 0, lpf_useIIR); %low-pass, 100Hz cutoff
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[filename_base '-lowpass'],'overwrite','on','gui','off'); 
    EEG = eeg_checkset( EEG );
else
    lpf = [];
end

%% Resample data to new_fs Hz
if new_fs < raw_fs
    EEG = pop_resample_myfilter(EEG,new_fs,0); % do not apply filter
    % EEG = pop_resample_myfilter( EEG, new_fs, lpf); % apply the same filter we used before
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[filename_base '-resampled'],'overwrite','on','gui','off'); 
    EEG = eeg_checkset( EEG );
end

%% Filter the data further
if ~isnan(highpass_cutoff) && highpass_cutoff < lowpass_cutoff
    [EEG,~,hpf] = pop_eegfilt_returnfilter( EEG, highpass_cutoff, 0, [], [0]); %high-pass, 0.5Hz cutoff
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[filename_base '-highpass'],'overwrite','on','gui','off'); 
    EEG = eeg_checkset( EEG );
else
    hpf = [];
end
if notch_bounds(2)>notch_bounds(1)
    [EEG,~,nf] = pop_eegfilt_returnfilter( EEG, notch_bounds(1), notch_bounds(2), [], [1], [1]); %notch, 59-61Hz, FFT
else
    warning('DJ:NotchBounds','Skipping Notch Filter!')
    nf = [];
end
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[filename_base '-filtered'],'overwrite','on','gui','off'); 

%% Add filters to EEGLAB dataset
EEG.etc.lowpassfilter = lpf;
EEG.etc.highpassfilter = hpf;
EEG.etc.notchfilter = nf;
%% Save the data
EEG = pop_saveset( EEG, 'filename',final_filename, 'filepath',data_dir);
disp('Data has been saved!');

%% Plot the data
%pop_eegplot( EEG, 1, 1, 1); % view the processed data

%% Clean up
toc % Display elapsed time
eeglab redraw