function y = ApplyHybridHdcaClassifier(data,trainingwindowlength,trainingwindowoffset,level2data,w,v)

% y = ApplyHybridHdcaClassifier(data,trainingwindowlength,trainingwindowoffset,level2data,w,v)
%
% INPUTS:
% - data is a [DxTxN] matrix of data, where D is # channels, T is # time
%   samples, and N is # trials.
% - trainingwindowlength is a scalar indicating the number of samples that
%   should be in each time window.
% - trainingwindowoffset is an M-element vector indicating the offset of 
%   each time window in samples.
% - level2data is a NxQ matrix of non-time-dependent data for each trial.
%   That is, each column of level2data will receive a "temporal" weight as
%   if it were another output from the spatial weights. 
% - w is a DxM matrix, where D is the number of electrodes and M is the
%   number of windows. w(:,i,j) is the set of spatial weights found by the 
%   classifier that best discriminates the training data in window i.
% - v is a [1x(M+Q)] matrix, in which v is the set of temporal/cross-modal
%   weights found by the classifier that best discriminates the data in 
%   the training set.
%
% OUTPUTS:
% - y is an Nx1 vector of 'interest scores', where a higher value of y(i)
% indicates that trial i is more likely to be a target.
%
% Created 3/11/14 by DJ.


nWindows = numel(trainingwindowoffset);
nTrials = size(data,3);
y_level1 = nan(nTrials,nWindows);

% Apply spatial weights to level 1 data
for iWin = 1:nWindows
    isInWin = trainingwindowoffset(iWin)+(0:trainingwindowlength-1);
    windowAvg = permute(mean(data(:,isInWin,:),2), [3,1,2]);
    y_level1(:,iWin) = windowAvg * w(:,iWin);
end

% Apply temporal weights to level 1 and 2 data
y_all = cat(2,y_level1,level2data);
y = y_all*v';