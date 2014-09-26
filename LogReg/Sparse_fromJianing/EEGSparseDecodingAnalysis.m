%%% -------------------------------------------------------------------%%%
%%% Recipe for Biosemi EEG
%%%
%%% Function: SingleTrialAnalysis
%%% This function performs single trial analysis for EEG data.
%%%
%%% Input: nothing
%%%     
%%% Output: nothing
%%% 
%%% Jianing Shi
%%% 09/13/2010
%%%
%%% -------------------------------------------------------------------%%%

%eeglab, close;

%%% Set up parameters related to single-trial analysis
% stimulus types
eegb.stimtype  = {'targapp','distapp'};
% analysis type: 'stim' - stimulus locked  'resp' - response locked
eegb.analtype  = 'stim';
% regressors type: 'binary'
eegb.regstype  = 'binary';
% classification types
eegb.bisetlist = [1,2];
% label type
eegb.labeltype = 'stim';
%efmr.labeltype = 'choice';

%%% classification initializations
% perform leave-one-out cross validation
eegb.LOO = 1; 	
% plot performance vs time (or frequency)	
eegb.showplot = 0; 	
% domain name for classification
eegb.domainname = 'time'; % or 'frequency'

switch eegb.domainname,
    case 'time',
        eegb.unit = 'ms'; % for labeling plots
        eegb.windowlength = 25; % width of training window
        eegb.windowoffset = 0:25:600; % centers of training window relative to event
    case 'frequency',
        eegb.unit = 'Hz'; % for labeling plots
        eegb.windowlength = 3;	% width of training window
        eegb.windowoffset = 0:100; % centers of training window relative to event
end


% Load EEG data for different stimuli types
for sti = 1:length(eegb.stimtype)

    fname = ['8-',eegb.stimtype{sti},'.set'];
    EEG = pop_loadset( 'filename', fname, 'filepath', ['./']);

    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, sti);

end % sti

% --------- binary classification --------------
if strcmp(eegb.regstype,'binary')

    % Binary classification for all three cases
    for cls = 1:size(eegb.bisetlist,1)

        classlabel = [eegb.stimtype{eegb.bisetlist(cls,1)} ...
                      eegb.stimtype{eegb.bisetlist(cls,2)}];

        chanexclude = [];

        LR = SparseLogisticRegression(ALLEEG,eegb.bisetlist(cls,:),classlabel,eegb.domainname,eegb.windowlength,eegb.windowoffset,eegb.labeltype,chanexclude);

    end % cls

end % binary classification

figure;
imagesc(LR.windowoffset,LR.chansubset,reshape(LR.wloo,LR.numchans,LR.numwins));
xlabel('time (ms)');
ylabel('channel');

load('biosemi64chanlocs.mat');
figure;
for widx = 1:length(LR.windowoffset)
    numrow = ceil(sqrt(length(LR.windowoffset)));
    subplot(numrow,numrow,widx);
    topoplot(LR.a{widx},chanlocs);
    title([num2str(LR.windowoffset(widx)),' ms']);
end % widx