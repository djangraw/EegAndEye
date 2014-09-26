function [y, w, v, fwdModel] = run_rsvp_classifier(ALLEEG,trainingwindowlength,trainingwindowoffset,cvmode,useica)

% Run a 2-step classifier similar to that used in LIINC RSVP experiments.
%
% [y, w, v, fwdModel] = run_rsvp_classifier(ALLEEG,trainingwindowlength,
%    trainingwindowoffset,cvmode,useica)
%
% INPUTS:
% - ALLEEG is a 2-element vector of EEGLAB data structures.  The first
% should contain distractor trials, and the second should contain target
% trials.
% - trainingwindowlength is a scalar indicating the number of samples that
% should be in each training window.
% -trainingwindowoffset is a vector indicating the offset of each training
% window in samples.
% -cvmode is a string specifying what type of cross-validation to run.  Can be:
%		'nocrossval' - run full model (no cross-validation)
%		'loo' - run leave-one-out cross-validation
%		'XXfold', where XX is an integer, will run XX-fold cross-validation
%	(DEFAULT:  '10fold')
% -useica is a binary value indicating whether the ica activations should
% be used as inputs (by default, the channel voltages will be used).
%
% OUTPUTS:
% - y is a n-element vector, where n is the number of trials in ALLEEG(1) 
% and ALLEEG(2) combined.  Each element is the "interest score" for that
% trial, a higher number indicating that the trial is more likely a target.
% - w is a DxMxP matrix, where D is the number of electrodes, M is the
% number of windows, and P is the number of folds in cross-validation.
% w(:,i,j) is the set of spatial weights found by the FLD that best 
% discriminates the data in window i in the training trials of fold j.
% - v is a 1xMxP matrix, in which v(:,:,j) is the set of temporal weights
% found by LR that best discriminates the data in the training trials of
% fold j.
% - fwdModel is a DxMxP matrix, where D is the number of electrodes, M is 
% the number of windows, and P is the number of folds in cross-validation.
% fwdModel(:,i,j) is the forward model taken by multiplying the inverse of
% the y values from window i and fold j by all the data from window i.
%
% Created 11/11/11 by DJ.
% Updated 12/13/11 by DJ - fixed y output bug
% Updated 3/29/13 by DJ - added useica input/option
% Updated 5/15/13 by DJ - updated fwdModels to work with ICs.
% Updated 6/3/13 by DJ - comments.

if nargin<4
    cvmode = '10fold';
end
if nargin<5 || isempty(useica)
    useica = false;
end

% Set up
cv = setGroupedCrossValidationStruct(cvmode,ALLEEG(1),ALLEEG(2));
nFolds = cv.numFolds;
nWindows = numel(trainingwindowoffset);

if useica
    nElecs = size(ALLEEG(1).icaact,1);
else
    nElecs = ALLEEG(1).nbchan;
end
nSamples = ALLEEG(1).pnts;
nTrials = ALLEEG(1).trials + ALLEEG(2).trials;
% Initialize
w = nan(nElecs,nWindows,nFolds);
v = nan(1,nWindows,nFolds);
y = nan(1,nTrials);
truth_trials = [zeros(ALLEEG(1).trials,1); ones(ALLEEG(2).trials,1)];
% Initialize weights y values
yTrain = cell(1,nFolds);
yTest = cell(1,nFolds);    
testTruthTrials = cell(1,nFolds);

for foldNum=1:nFolds   
    fprintf('Fold %d...\n',foldNum);
    % Separate out training and testing data
    if useica
        foldTrainingData = cat(3,ALLEEG(1).icaact(:,:,cv.incTrials1{foldNum}), ALLEEG(2).icaact(:,:,cv.incTrials2{foldNum}));
        foldTestingData = cat(3,ALLEEG(1).icaact(:,:,cv.valTrials1{foldNum}), ALLEEG(2).icaact(:,:,cv.valTrials2{foldNum}));
    else
        foldTrainingData = cat(3,ALLEEG(1).data(:,:,cv.incTrials1{foldNum}), ALLEEG(2).data(:,:,cv.incTrials2{foldNum}));
        foldTestingData = cat(3,ALLEEG(1).data(:,:,cv.valTrials1{foldNum}), ALLEEG(2).data(:,:,cv.valTrials2{foldNum}));
    end
    % Make corresponding truth matrices
    foldTrainingTruth = cat(3,zeros(1,nSamples,length(cv.incTrials1{foldNum})), ones(1,nSamples,length(cv.incTrials2{foldNum}))); 
    foldTestingTruth= cat(3,zeros(1,nSamples,length(cv.valTrials1{foldNum})), ones(1,nSamples,length(cv.valTrials2{foldNum}))); 
    testTruthTrials{foldNum} = [zeros(length(cv.valTrials1{foldNum}),1); ones(length(cv.valTrials2{foldNum}),1)];
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
    % Create ALLEEG struct with y as data
    foldTrainingTruth_trials = foldTrainingTruth(1,1,:);
    y0data = permute(yTrain{foldNum}(foldTrainingTruth_trials==0,:),[2,3,1]);
    y1data = permute(yTrain{foldNum}(foldTrainingTruth_trials==1,:),[2,3,1]);
    assignin('base','y0data',y0data);
    assignin('base','y1data',y1data);
    foldALLEEG(1) = pop_importdata('data','y0data');
    foldALLEEG(1) = eeg_checkset(foldALLEEG(1));
    foldALLEEG(2) = pop_importdata('data','y1data');
    foldALLEEG(2) = eeg_checkset(foldALLEEG(2));
    % Run logistic regression to get temporal weights
    regularize=1;
    lambda=1e1;
    lambdasearch=0;
    eigvalratio=1e-4;
    vinit=zeros(nWindows+1,1);
    show = 0;
    LOO = 0;
    ALLEEGout = pop_logisticregression_fast(foldALLEEG,[1 2],1:nWindows,1:nWindows,...
        1,1,regularize,lambda,lambdasearch,eigvalratio,vinit,show,LOO);
    % Extract temporal weights and find YIS interest scores
    v(:,:,foldNum) = ALLEEGout(1).icaweights;
    YisTrain = v(1,1:nWindows,foldNum)*yTrain{foldNum}';
    YisTest = v(1,1:nWindows,foldNum)*yTest{foldNum}';
    
    y([cv.valTrials1{foldNum}, ALLEEG(1).trials+cv.valTrials2{foldNum}]) = YisTest;
    
end

% Print Az
Az = rocarea(y,truth_trials);
fprintf('%s Az value: %.3f\n',cvmode,Az);

% Get forward models
fwdModel = zeros(ALLEEG(1).nbchan,nWindows,nFolds); % don't use nElecs - that may refer to #IC's.
for foldNum=1:nFolds
    for iWin = 1:nWindows
        % Extract relevant data
        isInWin = trainingwindowoffset(iWin)+(0:trainingwindowlength-1);
        if useica
            trainData = cat(3,ALLEEG(1).icaact(:,isInWin,cv.valTrials1{foldNum}), ALLEEG(2).icaact(:,isInWin,cv.valTrials2{foldNum}));
            rawData = cat(3,ALLEEG(1).data(:,isInWin,cv.valTrials1{foldNum}), ALLEEG(2).data(:,isInWin,cv.valTrials2{foldNum}));
        else
            trainData = cat(3,ALLEEG(1).data(:,isInWin,cv.valTrials1{foldNum}), ALLEEG(2).data(:,isInWin,cv.valTrials2{foldNum}));
            rawData = trainData;
        end
        % Get y values
        dataAvg = permute(mean(trainData,2), [3,1,2]);
        yAvg = dataAvg*w(:,iWin,foldNum);
        % Compute forward model        
        rawDataAvg = permute(mean(rawData,2), [3,1,2]);
        fwdModel(:,iWin,foldNum) = yAvg \ rawDataAvg;
    end
end

% % PLOTTING: Scatter plot
% y_all = cat(1,yTest{:});
% ist_all = cat(1,testTruthTrials{:})==1;
% scatter(y_all(~ist_all,1),y_all(~ist_all,2),'b.')
% hold on
% scatter(y_all(ist_all,1),y_all(ist_all,2),'r.')
% legend('distractors','targets')
% xlabel(sprintf('y(%g)',ALLEEG(1).times(trainingwindowoffset(1))))
% ylabel(sprintf('y(%g)',ALLEEG(1).times(trainingwindowoffset(2))))
% title(sprintf('%s Az value: %.3f',cvmode,Az));
