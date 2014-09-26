function [y, w, v, fwdModel, y_level1] = TrainHybridHdcaClassifier(data, truth, trainingwindowlength, trainingwindowoffset, level2data, fwdModelData,cvmode)

% Run a hierarchical classifier with added "level 2" data.
%
% [y, w, v, fwdModel, y_level1] = TrainHybridHdcaClassifier(data, truth, trainingwindowlength,
%    trainingwindowoffset, level2data, fwdModelData)
%
% INPUTS:
% - data is a [DxTxN] matrix of data, where D is # channels, T is # time,
%   and N is # trials.
% - truth is a 1xN matrix of binary labels indicating the class of each
%   trial.
% - trainingwindowlength is a scalar indicating the number of samples that
%   should be in each training window.
% - trainingwindowoffset is an M-element vector indicating the offset of 
%   each training window in samples.
% - level2data is a NxQ matrix of non-time-dependent data for each trial.
%   That is, each column of level2data will receive a "temporal" weight as
%   if it were another output from the spatial weights.
% - fwdModelData is a [ExTxN] matrix (where E is # channels) used to
%   calculate the forward models. (e.g., if ICA activations are 'data', the
%   corresponding raw data could be used as 'fwdModelData').
% - cvmode is a string used as input to setCrossValidationStruct.
%  'nocrossval' means all data is used during training. [default: '10fold']
%
% OUTPUTS:
% - y is a N-element vector indicating the "interest score" for each trial,
%   where a higher number indicates that the trial is more likely a target.
% - w is a DxM matrix, where D is the number of electrodes and M is the
%   number of windows.
%   w(:,i) is the set of spatial weights found by the FLD that best 
%   discriminates the data in window i in the training trials.
% - v is a [1x(M+Q)] matrix, in which v is the set of temporal weights 
%   found by LR that best discriminates the data in the training trials.
% - fwdModel is a DxM matrix, where D is the number of electrodes and M is 
%   the number of windows.
%   fwdModel(:,i) is the forward model taken by multiplying the inverse 
%   of the y values from window i by all the data from window i.
% - y_level1 is an NxM matrix of the 'within-bin interest scores' for each 
%   bin of the EEG.
%
% Created 3/11/14 by DJ based on RunHybridHdcaClassifier.
% Updated 9/25/14 by DJ - empty fwdModelData defaults to data.

% Handle defaults
if ~exist('level2data','var')
    level2data = [];
end
if ~exist('fwdModelData','var') || isempty(fwdModelData)
    fwdModelData = data;
end
if ~exist('cvmode','var') || isempty(cvmode)
    cvmode = '10fold'; % for getting evaluation set trials
end

% Set up
nWindows = numel(trainingwindowoffset);
nAddOns = size(level2data,2);
nWindowWeights = nWindows + nAddOns;
[nElecs, nSamples, nTrials] = size(data);
nFmElecs = size(fwdModelData,1); % for calculating forward model

% Make sample-by-sample truth matrix
sampleTruth = repmat(reshape(truth,[1 1 nTrials]),[1 nSamples 1]);

% Set up logistic regression params
params.regularize=1;
params.lambda=1e1;
params.lambdasearch=0;
params.eigvalratio=1e-4;
params.vinit=zeros(nWindowWeights+1,1);
params.show = 0;
params.LOO = 0;    

%%% TRAIN LEVEL 1 %%%
tic;

cv_fold = setCrossValidationStruct(cvmode,nTrials);
wFold = zeros(nElecs,nWindows,cv_fold.numFolds);
yEval = zeros(nTrials,nWindows);

% INNER LOOP
for jFold = 1:cv_fold.numFolds
    innerTrainingData = data(:,:,cv_fold.incTrials{jFold});
    innerTestingData = data(:,:,cv_fold.valTrials{jFold});
    innerTrainingTruth = sampleTruth(:,:,cv_fold.incTrials{jFold});
    for iWin=1:nWindows
        % Extract relevant data
        isInWin = trainingwindowoffset(iWin)+(0:trainingwindowlength-1);
        trainingData = permute(reshape(innerTrainingData(:,isInWin,:),size(innerTrainingData,1),length(isInWin)*size(innerTrainingData,3)), [2 1]);    
        testingData = permute(reshape(innerTestingData(:,isInWin,:),size(innerTestingData,1),length(isInWin)*size(innerTestingData,3)), [2 1]);    
        trainingTruth = permute(reshape(innerTrainingTruth(:,isInWin,:),size(innerTrainingTruth,1),length(isInWin)*size(innerTrainingTruth,3)), [2 1]);
        % Find spatial weights using FLD
        try
            [~,~,~,~,coeff] = classify(testingData,trainingData,trainingTruth);
        catch
            warning('classify errored... trying diaglinear option.')
            [~,~,~,~,coeff] = classify(testingData,trainingData,trainingTruth,'diaglinear');
        end
        wFold(:,iWin,jFold) = coeff(2).linear;
%             w(:,iWin,foldNum,jFold) = coeff(2).linear;
        % Get 'evaluation' y values
        testingAvg = permute(mean(innerTestingData(:,isInWin,:),2), [3,1,2]);         
        yEval(cv_fold.valTrials{jFold},iWin) = testingAvg*wFold(:,iWin,jFold);            
    end        
end
% Get mean weights across folds
w = mean(wFold,3);
yTest = nan(nTrials,nWindows);
% use avg of wts across inner folds to find 'testing y values'
for iWin = 1:nWindows
    isInWin = trainingwindowoffset(iWin)+(0:trainingwindowlength-1);
    testingAvg = permute(mean(data(:,isInWin,:),2), [3,1,2]);
    yTest(:,iWin) = testingAvg * w(:,iWin);
end


%%% TRAIN LEVEL 2 %%%

% Get add-on (level 2) data
yEval_AddOn = zeros(nTrials,nAddOns);    
wAddOn = nan(1,nAddOns);
for iAddOn = 1:nAddOns
    if all(level2data(:,iAddOn)==1) % if this column is a simple offset...
        yEval_AddOn(:,iAddOn) = 1;
    else
        testingData_AddOn = level2data(:,iAddOn);
        % Perform dummy classification to scale data properly
        [~,~,~,~,coeff] = classify(testingData_AddOn, trainingData_AddOn,foldTrainingTruth_trials(:));
        wAddOn(iAddOn) = coeff(2).linear;
        if wAddOn(iAddOn)==0 % Added for special case where class means are equal
            wAddOn(iAddOn) = eps;
        end
        yEval_AddOn(:,iAddOn) = testingData_AddOn*wAddOn(iAddOn);
    end
end
yTest_AddOn = yEval_AddOn;

% Append level2data to yEval and yTest
yEvalData = cat(2,yEval,yEval_AddOn);
yTestData = cat(2,yTest,yTest_AddOn);

% Standardize stddev of each bin
stdEval = std(yEvalData,0,1);
for iWin=1:nWindowWeights        
    yTestData(:,iWin) = yTestData(:,iWin)/stdEval(iWin); % normalize using eval data
    yEvalData(:,iWin) = yEvalData(:,iWin)/stdEval(iWin);
end

% Run logistic regression to get temporal weights
[~,~,stats] = RunSingleLR(permute(yEvalData,[2,3,1]),truth,params);
% Extract temporal weights and find YIS interest scores
v0 = stats.wts(1:nWindowWeights)';
y = v0 * yTestData';

% Assemble v vector to incorporate FLD and std corrections
v = v0;
v((nWindows+1):end) = v((nWindows+1):end) .* wAddOn;
v = v ./ stdEval;

% Report elapsed time
tRun = toc;
fprintf('Done with training! Took %.1f seconds.\n',tRun);
    

%%% CLEAN UP %%%

% set level1 y values
y_level1 = yTest;


% Print Az
Az = rocarea(y,truth);
fprintf('Training Az value: %.3f\n', Az);

% Get forward models
fwdModel = zeros(nFmElecs,nWindows);
for iWin = 1:nWindows
    % Extract relevant data
    isInWin = trainingwindowoffset(iWin)+(0:trainingwindowlength-1);
    testData = data(:,isInWin,:);
    fmData = fwdModelData(:,isInWin,:);
    % Get y values
    dataAvg = permute(mean(testData,2), [3,1,2]);
    fmDataAvg = permute(mean(fmData,2), [3,1,2]);
    yAvg = dataAvg * w(:,iWin);
    % Compute forward model
    fwdModel(:,iWin) = yAvg \ fmDataAvg;
end

