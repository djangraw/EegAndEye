function PlotSingleTrials(EEG,channels)

% PlotSingleTrials(EEG,channels)
%
% INPUTS:
% -EEG is an eeglab data structure of the same name.
% -channels is a vector of the channel numbers you want to see, or a vector
% of cells containing the names of channels you want to see.
%
% Created 12/1/10 by DJ.


%% Convert to numbers, if necessary
chanLabels = {EEG.chanlocs(:).labels}; % the names of all the channels
if iscell(channels)
    channels = find(ismember(chanLabels,channels));
end

clear chanLabels

%% Plot each 
for i=1:numel(channels)
    pop_prop( EEG, 1, channels(i), NaN, {'freqrange' [2 50] });
end