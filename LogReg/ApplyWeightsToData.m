function [y, xw] = ApplyWeightsToData(EEG,trainingwindowlength,trainingwindowoffset,w,v);

% Apply classifier weights to data
%
% [y, xw] =
% ApplyWeightsToData(EEG,trainingwindowlength,trainingwindowoffset,w,v);
%
% INPUTS:
% - EEG is an epoched EEGLAB dataset.
% - trainingwindowlength is a scalar indicating the number of samples that
% should be in each training window.
% -trainingwindowoffset is an M-element vector indicating the offset of 
% each training window in samples.
% - w is a matrix of the spatial weights output by a 2-level classifier.
% w is a DxMxP matrix, where D is the number of electrodes, M is the
% number of windows, and P is the number of folds in cross-validation.
% w(:,i,j) is the set of spatial weights found by the FLD that best 
% discriminates the data in window i in the training trials of fold j.
% - v is a matrix of the temporal weights output by a 2-level classifier.
% v is a [1x(M+Q)xP] matrix, in which v(:,:,j) is the set of temporal weights
% found by LR that best discriminates the data in the training trials of
% fold j.
%
% OUTPUTS:
% - y is an N-element vector (where N is the number of trials) of the y 
% values resulting from applying w and v to the testing data. 
% - xw is an NxM matrix of the data in each bin multiplied by the spatial
% weights (w) for that bin.
%
% Created 10/30/13 by DJ.

% Set up
w_mean = mean(w,3); % elecs x bins 
v_mean = mean(v,3); % 1 x bins
nBins = length(trainingwindowoffset);
xw = nan(EEG.trials,nBins);

% Main loop
for i=1:nBins
    fprintf('Bin %d...\n',i);
    % Find relevant data
    iInWin = trainingwindowoffset(i)-1+(1:trainingwindowlength);
    x = permute(mean(EEG.data(:,iInWin,:),2),[3 1 2]); % trials x elecs x 1
    % Apply spatial weights
    xw(:,i) = x * w_mean(:,i); % trials x 1
end    
% Apply temporal weights
y = xw * v_mean'; 

