function EEG = NEDE_FilterEegData(eegdata, prefix,subject,session, raw_fs,new_fs, output_suffix,bandpass_bounds,notch_bounds, electrodeloc_filename)

% EEG = NEDE_FilterEegData(eegdata, prefix,subject,session, raw_fs,new_fs, output_suffix,bandpass_bounds,notch_bounds, electrodeloc_filename)
%
% Created 9/25/14 by DJ based on NEDE_ImportEegData.

if ~exist('electrodeloc_filename','var') || isempty(electrodeloc_filename)
    if size(eegdata,1)==9
        electrodeloc_filename = 'dave_abm9chan.loc';
    end
end


if isempty(new_fs)
    new_fs = old_fs;
end
% Set LPF preferences
lpf_useIIR = 0;
lpf_order = [];
% Set filter bounds
highpass_cutoff = bandpass_bounds(1);
lowpass_cutoff = bandpass_bounds(2);

% Directories
data_dir = [cd '/'];
electrodeloc_dir = '/Users/dave/Documents/Tools/eeglab/locationfiles/';

% Filenames
filename_base = sprintf('%s-%d-%d',prefix,subject,session);
final_filename = sprintf('%s%s.set',filename_base,output_suffix); 


% Import Data
EEG = pop_importdata('dataformat','array','nbchan',0,'srate',raw_fs,...
    'subject',num2str(subject),'pnts',0,'xmin',0,'session',num2str(session));
EEG.data = eegdata;

%
EEG.setname = sprintf('%s-%d-%d-raw',prefix,subject,session);
EEG = eeg_checkset( EEG );

%% Put in the channel location information
EEG=pop_chanedit(EEG, 'load',{[electrodeloc_dir electrodeloc_filename] 'filetype' 'autodetect'});
EEG = eeg_checkset( EEG );

%% Low-Pass Filter the data
if ~isnan(lowpass_cutoff) && highpass_cutoff < lowpass_cutoff
    [EEG, ~, lpf] = pop_eegfilt_returnfilter( EEG, 0, lowpass_cutoff, lpf_order, 0, lpf_useIIR); %low-pass, 100Hz cutoff
    EEG = eeg_checkset( EEG );
else
    lpf = [];
end

%% Resample data to new_fs Hz
if new_fs < raw_fs
    EEG = pop_resample_myfilter(EEG,new_fs,0); % do not apply filter
    EEG = eeg_checkset( EEG );
end

%% Filter the data further
if ~isnan(highpass_cutoff) && highpass_cutoff < lowpass_cutoff
    [EEG,~,hpf] = pop_eegfilt_returnfilter( EEG, highpass_cutoff, 0, [], 0); %high-pass, 0.5Hz cutoff
    EEG = eeg_checkset( EEG );
else
    hpf = [];
end
if notch_bounds(2)>notch_bounds(1)
    [EEG,~,nf] = pop_eegfilt_returnfilter( EEG, notch_bounds(1), notch_bounds(2), [], 1, 1); %notch, 59-61Hz, FFT
else
    warning('DJ:NotchBounds','Skipping Notch Filter!')
    nf = [];
end

%% Add filters to EEGLAB dataset
EEG.etc.lowpassfilter = lpf;
EEG.etc.highpassfilter = hpf;
EEG.etc.notchfilter = nf;
%% Save the data
EEG.setname = sprintf('%s%s',filename_base,output_suffix); % rename set
EEG = pop_saveset( EEG, 'filename',final_filename, 'filepath',data_dir);
disp('Data has been saved!');

%% Plot the data
%pop_eegplot( EEG, 1, 1, 1); % view the processed data

%% Clean up
toc % Display elapsed time
% eeglab redraw