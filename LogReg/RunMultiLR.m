function [Az, AzLoo, stats] = RunMultiLR(data, trialTruth, windowlength, windowslide)

% Run logistic regression and leave-one-out analysis on general data.
%
% [Az, AzLoo] = RunMultiLR(data, trialTruth, windowlength, windowslide)
%
% INPUTS:
% - data is an dxmxn matrix where each row is a feature, each column is an
% IID sample from that feature, and each page is a trial.
% - trialTruth is an n-element vector containing binary labels for the class of
% each sample.
% - windowlength is the number of samples long that you want your window to
% be.
% - windowslide is the number of samples you want your window to slide each
% time.
%
% OUTPUTS:
% - Az is the area under the ROC curve for the training data for each window.
% - AzLoo is the LOO Az value for each window.
% - stats is a struct for each window containing fields wts, fwdModel, y, 
% wtsLoo, fwdModelLoo, and yLoo.
%
% Created 1/20/11 by DJ.
% Updated 9/11/12 by DJ - added stats output

if nargin<3 || isempty(windowlength)
    windowlength = 13;
end
if nargin<4 || isempty(windowslide)
    windowslide = 3;
end

nSamples = size(data,2);
nSlides = floor((nSamples-windowlength)/windowslide);

Az = NaN(1,nSlides);
AzLoo = NaN(1,nSlides);
stats = cell(1,nSlides);
parfor i=1:nSlides
    fprintf('===Running window %d of %d\n',i,nSlides);
    thisSamples = (i-1)*windowslide + (1:windowlength);
    thisData = data(:,thisSamples,:);
    [Az(i), AzLoo(i), stats{i}] = RunSingleLR(thisData,trialTruth);
end
% Turn cell array into struct vector
try
    stats = [stats{:}];
end
