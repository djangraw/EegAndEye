function [y, w, v, fwdModel] = run_rsvp_classifier_rawdata(data, truth, trainingwindowlength, trainingwindowoffset, cvmode, level2data, fwdModelData)

% Run a 2-step classifier similar to that used in LIINC RSVP experiments.
%
% [y, w, v, fwdModel] = run_rsvp_classifier_rawdata(data, truth, trainingwindowlength,
%    trainingwindowoffset,cvmode,level2data)
%
% INPUTS:
% - data is a [DxTxN] matrix of data, where D is # channels, T is # time,
% and N is # trials.
% - truth is a 1xN matrix of binary labels indicating the class of each
% trial.
% - trainingwindowlength is a scalar indicating the number of samples that
% should be in each training window.
% -trainingwindowoffset is an M-element vector indicating the offset of 
% each training window in samples.
% -cvmode is a string indicating the cross-validation mode, used as input
% to setGroupedCrossValidationStruct.m.
% -level2data is a NxQ matrix of non-time-dependent data for each trial.
% That is, each column of level2data will receive a "temporal" weight as if
% it were another output from the spatial weights.
% -fwdModelData is a [ExTxN] matrix (where E is # channels) used to
% calculate the forward models.
%
% OUTPUTS:
% - y is a N-element vector, indicating the "interest score" for each
% trial, a higher number indicating that the trial is more likely a target.
% - w is a DxMxP matrix, where D is the number of electrodes, M is the
% number of windows, and P is the number of folds in cross-validation.
% w(:,i,j) is the set of spatial weights found by the FLD that best 
% discriminates the data in window i in the training trials of fold j.
% - v is a [1x(M+Q)xP] matrix, in which v(:,:,j) is the set of temporal weights
% found by LR that best discriminates the data in the training trials of
% fold j.
% - fwdModel is a DxMxP matrix, where D is the number of electrodes, M is 
% the number of windows, and P is the number of folds in cross-validation.
% fwdModel(:,i,j) is the forward model taken by multiplying the inverse of
% the y values from window i and fold j by all the data from window i.
%
% Created 11/11/11 by DJ.
% Updated 12/13/11 by DJ - fixed y output bug
% Updated 3/25/13 by DJ - _rawdata, cv input

% Handle defaults
if ~exist('level2data','var')
    level2data = [];
end
if ~exist('fwdModelData','var');
    fwdModelData = data;
end

% Set up
nWindows = numel(trainingwindowoffset);
nWindowWeights = nWindows + size(level2data,2);
[nElecs, nSamples, nTrials] = size(data);
nFmElecs = size(fwdModelData,1); % for calculating forward model
cv = setCrossValidationStruct(cvmode,nTrials);
nFolds = cv.numFolds;
% Initialize
w = nan(nElecs,nWindows,nFolds);
v = nan(1,nWindowWeights,nFolds);
y = nan(1,nTrials);
% Initialize weights y values
yTrain = cell(1,nFolds);
yTest = cell(1,nFolds);    
sampleTruth = repmat(reshape(truth,[1 1 nTrials]),[1 nSamples 1]);
testTruth = cell(1,nFolds);

for foldNum=1:nFolds   
    fprintf('Fold %d...\n',foldNum);
    % Separate out training and testing data
    foldTrainingData = data(:,:,cv.incTrials{foldNum});
    foldTestingData = data(:,:,cv.valTrials{foldNum});
    % Make corresponding truth matrices
    foldTrainingTruth = sampleTruth(:,:,cv.incTrials{foldNum});    
    testTruth{foldNum} = truth(cv.valTrials{foldNum});
    foldTrainingTruth_trials = foldTrainingTruth(1,1,:); % used later
    % Initialize weights y values
    yTrain{foldNum} = nan(size(foldTrainingData,3),nWindows);
    yTest{foldNum} = nan(size(foldTestingData,3),nWindows);  
    for iWin=1:nWindows
        % Extract relevant data
        isInWin = trainingwindowoffset(iWin)+(0:trainingwindowlength-1);
        trainingData = permute(reshape(foldTrainingData(:,isInWin,:),size(foldTrainingData,1),length(isInWin)*size(foldTrainingData,3)), [2 1]);    
        testingData = permute(reshape(foldTestingData(:,isInWin,:),size(foldTestingData,1),length(isInWin)*size(foldTestingData,3)), [2 1]);    
        trainingTruth = permute(reshape(foldTrainingTruth(:,isInWin,:),size(foldTrainingTruth,1),length(isInWin)*size(foldTrainingTruth,3)), [2 1]);
        % Find spatial weights using FLD
        [~,~,~,~,coeff] = classify(testingData,trainingData,trainingTruth);
        w(:,iWin,foldNum) = coeff(2).linear;
        % Get y values
        trainingAvg = permute(mean(foldTrainingData(:,isInWin,:),2), [3,1,2]);
        testingAvg = permute(mean(foldTestingData(:,isInWin,:),2), [3,1,2]);
        yTrain{foldNum}(:,iWin) = trainingAvg*w(:,iWin,foldNum);
        yTest{foldNum}(:,iWin) = testingAvg*w(:,iWin,foldNum);
    end
    
    % Get add-on data
    yTrain_AddOn = zeros(length(cv.incTrials{foldNum}),size(level2data,2));
    yTest_AddOn = zeros(length(cv.valTrials{foldNum}),size(level2data,2));
    for iAddOn = 1:size(level2data,2)
        if all(level2data(:,iAddOn)==1) % for an offset
            yTrain_AddOn(:,iAddOn) = 1;
            yTest_AddOn(:,iAddOn) = 1;
        else
            trainingData_AddOn = level2data(cv.incTrials{foldNum},iAddOn);
            testingData_AddOn = level2data(cv.valTrials{foldNum},iAddOn);
            % Perform dummy classification to scale data properly
            [~,~,~,~,coeff] = classify(testingData_AddOn, trainingData_AddOn,foldTrainingTruth_trials(:));
            wAddOn = coeff(2).linear;
            yTrain_AddOn(:,iAddOn) = trainingData_AddOn*wAddOn;
            yTest_AddOn(:,iAddOn) = testingData_AddOn*wAddOn;
        end
    end
    % Append level2data to yTrain and yTest
    yTrainData = cat(2,yTrain{foldNum},yTrain_AddOn);
    yTestData = cat(2,yTest{foldNum},yTest_AddOn);
    
    % Standardize stddev of each bin
    for iWin=1:nWindowWeights
        yTestData(:,iWin) = yTestData(:,iWin)/std(yTrainData(:,iWin)); % normalize using training data
        yTrainData(:,iWin) = yTrainData(:,iWin)/std(yTrainData(:,iWin));        
    end
    
    % Create ALLEEG struct with y as data
    y0data = permute(yTrainData(foldTrainingTruth_trials==0,:),[2,3,1]);
    y1data = permute(yTrainData(foldTrainingTruth_trials==1,:),[2,3,1]);
    assignin('base','y0data',y0data);
    assignin('base','y1data',y1data);
    foldALLEEG(1) = pop_importdata('data','y0data');
    foldALLEEG(1) = eeg_checkset(foldALLEEG(1));
    foldALLEEG(2) = pop_importdata('data','y1data');
    foldALLEEG(2) = eeg_checkset(foldALLEEG(2));
    % Run logistic regression to get temporal weights
    regularize=1;
%     lambda = [repmat(nElecs,1,nWindows), ones(1,nWindowWeights-nWindows)];
    lambda=1e1;
    lambdasearch=0;
    eigvalratio=1e-4;
    vinit=zeros(nWindowWeights+1,1);
    show = 0;
    LOO = 0;
    ALLEEGout = pop_logisticregression_fast(foldALLEEG,[1 2],1:nWindowWeights,1:nWindowWeights,...
        1,1,regularize,lambda,lambdasearch,eigvalratio,vinit,show,LOO);
    % Extract temporal weights and find YIS interest scores
    v(:,:,foldNum) = ALLEEGout(1).icaweights;
    YisTrain = v(1,1:nWindowWeights,foldNum)*yTrainData';
    YisTest = v(1,1:nWindowWeights,foldNum)*yTestData';
    
    y(cv.valTrials{foldNum}) = YisTest;
    
end

% Print Az
Az = rocarea(y,truth);
if cv.numFolds==nTrials
    cvmode = 'loo';
else
    cvmode = sprintf('%d-fold',cv.numFolds);
end
fprintf('%s Az value: %.3f\n',cvmode,Az);

% Get forward models
fwdModel = zeros(nFmElecs,nWindows,nFolds);
for foldNum=1:nFolds
    for iWin = 1:nWindows
        % Extract relevant data
        isInWin = trainingwindowoffset(iWin)+(0:trainingwindowlength-1);
        trainData = data(:,isInWin,cv.valTrials{foldNum});
        fmData = fwdModelData(:,isInWin,cv.valTrials{foldNum});
        % Get y values
        dataAvg = permute(mean(trainData,2), [3,1,2]);
        fmDataAvg = permute(mean(fmData,2), [3,1,2]);
        yAvg = dataAvg*w(:,iWin,foldNum);
        % Compute forward model
        fwdModel(:,iWin,foldNum) = yAvg \ fmDataAvg;
    end
end

% % TEMP: Scatter plot
% y_all = cat(1,yTest{:});
% ist_all = cat(1,testTruth{:})==1;
% scatter(y_all(~ist_all,1),y_all(~ist_all,2),'b.')
% hold on
% scatter(y_all(ist_all,1),y_all(ist_all,2),'r.')
% legend('distractors','targets')
% xlabel(sprintf('y(%g)',ALLEEG(1).times(trainingwindowoffset(1))))
% ylabel(sprintf('y(%g)',ALLEEG(1).times(trainingwindowoffset(2))))
% title(sprintf('%s Az value: %.3f',cvmode,Az));
