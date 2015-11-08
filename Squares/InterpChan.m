function EEG = InterpChan(EEG, chanList, iTrials)

% EEG = InterpChan(EEG, chanList, iTrials)
%
% e.g., EEG = InterpChan(EEG, [4, 7, 15])
%       EEG = InterpChan(EEG, {'NZ','F9','F10'}, 25);
%       EEG = InterpChan(EEG, EEG.chanlocs([4, 7, 15]), 1:10);
%
% Given an EEG dataset, this function replaces the channels in chansList by
% interpolating it with all the other channels in the array. The affect of
% each channel is scaled by a factor of 1/distance from channel being
% interpolated
%
% Inputs
%   1. EEG: EEGLAB structure of the second dataset
%   2. chanList: vector of indices, labels, or chanlocs structs of channels
%   to be interpolated (see examples above)
%   3. iTrials: indices of trials on which data on the given electrodes 
%   should be interpolated [default: all trials]
%
% Outputs
%   1. EEG: modified EEG structure with interpolated electrodes
%
% Created 7/24/12 by Ansh Johri
% Modified 7/24/12 by Ansh Johri
% Modified 8/7/12 by DJ - comments, chanNumVec calculation, added iTrials
%  input
% Modified 8/8/12 by DJ - allow chanList to be list of labels
% Updated 2/25/14 by DJ - document which chans were interpolated

if nargin<3 || isempty(iTrials)
    iTrials = 1:EEG.trials; % default is all trials
end

% Go through each channel to be interpolated, given in input list 
n = numel(chanList);

% Convert input to vector of channel indices
if isstruct(chanList)
    [~, chanNumVec] = ismember({chanList.labels},{EEG.chanlocs.labels});
elseif iscell(chanList)
    [~, chanNumVec] = ismember(chanList,{EEG.chanlocs.labels});
elseif isnumeric(chanList)
    chanNumVec = chanList;
end

% Create a number vector containing the channel indices being interpolated

for i = 1:n
            
    % Find the coordinates of the current channel being interpolated.
    badX = EEG.chanlocs(chanNumVec(i)).X;
    badY = EEG.chanlocs(chanNumVec(i)).Y;
    badZ = EEG.chanlocs(chanNumVec(i)).Z; 
        
    % Go through every channel
    smallSum = 0;
    largeSum = 0;
    for j = 1: EEG.nbchan

        % Consider electrode only if it isn't one of the ones being
        % interpolated
        if ~any(chanNumVec == j)
            
            % Calculate angle along unit sphere between the two channels
            % using dot product method
            goodX = EEG.chanlocs(j).X;
            goodY = EEG.chanlocs(j).Y;
            goodZ = EEG.chanlocs(j).Z;

            dotProd = badX*goodX + badY*goodY + badZ*goodZ;
            theta = acos(dotProd) ;   
            
            % Calculate distance based on theta in the unit circle
            dist = theta;
            invDist = 1/dist;
            invDistAll(j) = invDist;
            % Weight this electrode's data by invDist
            newData = invDist * EEG.data(j,:,iTrials);
            
            % Compile sums
            smallSum = smallSum + invDist; % sum of distances
            largeSum = largeSum + newData; % sum of data

        end      
    end
    invDistAll = invDistAll/sum(invDistAll);
    % Determine final data to be interpolated into channel of interest
    EEG.data(chanNumVec(i),:,iTrials) = largeSum / smallSum;
end

% document which chans were interpolated
chanLabels = {EEG.chanlocs.labels};
EEG.etc.chansInterpolated = chanLabels(chanNumVec); 





 