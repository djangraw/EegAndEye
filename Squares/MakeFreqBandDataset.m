function EEG = MakeFreqBandDataset(subject,freqband,bandname,saveresult)

% Created 1/30/13 by DJ.

% declare defaults
if nargin<2 || isempty(freqband)
    freqband = [9 14];
end
if nargin<3 || isempty(bandname)
    bandname = 'alpha';
end
if nargin<4 || isempty(saveresult)
    saveresult = false;
end

% Load dataset
disp('Loading...')
EEG = pop_loadset('filename',sprintf('sq-%d-all-filtered-50Hz-interpduds.set',subject));
% Filter dataset
disp('Filtering...')
EEG = pop_iirfilt( EEG, freqband(1), freqband(2), [], 0);
% Take Hilbert transform
disp('Taking Hilbert Transform...')
hilb_data = hilbert(EEG.data')';
EEG.data = abs(hilb_data);
EEG = pop_editset(EEG, 'setname', sprintf('%s %s-IIR', EEG.setname, bandname));
% Save changes to dataset
if saveresult
    disp('Saving new dataset...')    
    EEG = pop_saveset(EEG,'filename',sprintf('sq-%d-all-filtered-50Hz-interpduds-%s.set',subject,bandname));
else
    disp('Returning new dataset (without saving)...');   
end
disp('Done!')