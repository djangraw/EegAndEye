function EEG = ImportAlphaDataset(file_in, file_out, remove_elecs,interp_elecs)

% EEG = ImportAlphaDataset(file_in, file_out, remove_elecs,interp_elecs)
%
% Created 3/18/13 by DJ.

% set up
scaling = 1e6;
% Import & save
eegdata = readEEG_sensorium2013(file_in);
eegdata = eegdata*scaling;
assignin('base','eegdata',eegdata); % send to base workspace where eeglab will find it.
EEG = pop_importdata('dataformat','array','nbchan',0,'data','eegdata','setname',file_out,'srate',1000,'pnts',0,'xmin',0,'chanlocs','/Users/dave/Documents/Tools/eeglab/locationfiles/dave_sensorium90_84chan.loc');
pop_saveset(EEG,'filename',file_out);
% Remove electrodes
RemoveElectrodes(file_out, [file_out '-noduds'], remove_elecs);
% Interpolate electrodes
InterpolateElectrodes([file_out '-noduds'], [file_out '-interpduds'], interp_elecs);