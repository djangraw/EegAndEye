% LoadAllEpochs.m
%
% Variable 'subjects' must be defined.
%
% Created 1/18/11 by DJ.

% Clear EEGLAB
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];

% Load each subject in turn
for i = 1:numel(subjects)
    subject = subjects(i);
    LoadEpochs;
end
% Update EEGLAB GUI
eeglab redraw;