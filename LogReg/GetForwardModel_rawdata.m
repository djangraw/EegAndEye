function [fwdModel, cv] = GetForwardModel_rawdata(data,y,trainingwindowoffset,trainingwindowlength,cvmode)

% fwdModel =
% GetForwardModel_rawdata(data,y,trainingwindowoffset,trainingwindowlength)
%
% INPUTS:
% -data is a [DxTxN] matrix (D = # chans, T = # samples, N = # trials)
% -y is a  [NxP] matrix (P = # time windows)
% -trainingwindowoffset is a P-element vector
% -trainingwindowlength is a scalar
% -cvmode is a string indicating the cross-validation mode (used as input
% to setCrossValidationStruct.m).
%
% OUTPUTS:
% -fwdModel is a [DxPxF] matrix (F = # folds)
%
% Created 5/23/13 by DJ.

if ~exist('cvmode','var')
    cvmode = 'nocrossval';
end

nWindows = length(trainingwindowoffset);
[nChan, ~, nTrials] = size(data,3);
cv = setCrossValidationStruct(cvmode,nTrials);
nFolds = cv.numFolds;

% Get forward models
fwdModel = zeros(nChan,nWindows,nFolds); % don't use nElecs - that may refer to #IC's.
for foldNum=1:nFolds
    for iWin = 1:nWindows
        % Extract relevant data
        isInWin = trainingwindowoffset(iWin)+(0:trainingwindowlength-1);
        rawData = data(:,isInWin,cv.valTrials{foldNum});
        % Extract relevant y values
        yAvg = y(cv.valTrials{foldNum},iWin);
        % Compute forward model        
        rawDataAvg = permute(mean(rawData,2), [3,1,2]);
        fwdModel(:,iWin,foldNum) = yAvg \ rawDataAvg;
    end
end
