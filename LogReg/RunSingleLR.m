function [Az, AzLoo, stats, AzLto] = RunSingleLR(data, trialTruth, params)

% Run logistic regression and leave-one-out analysis on general data.
%
% [Az, AzLoo, stats, AzLto] = RunSingleLR(data, trialTruth, params)
%
% INPUTS:
% - data is an dxmxn matrix where each row is a feature, each column is an
% IID sample from that feature, and each page is a trial.
% - trialTruth is an n-element vector containing binary labels for the class of
% each sample.
% - params is a struct with fields 'regularize', 'lambda', 'lambdasearch', 
% 'eigvalratio','vinit','show', and 'LOO'.
%
% OUTPUTS:
% - Az is the area under the ROC curve for the training data.
% - AzLoo is the LOO Az value.
% - stats is a struct containing fields wts, fwdModel, y, wtsLoo,
% fwdModelLoo, and yLoo.
% - AzLto is the LTO Az value.
%
% Created 1/4/11 by DJ.
% Updated 1/20/11 by DJ - changed name from RunLR, which was already taken
% by Jen's program.
% Updated 9/11/12 by DJ - added stats struct, LOO fwd model
% Updated 5/29/13 by DJ - added params input
% Updated 8/24/13 by DJ - removed for loops for efficiency
% Updated 12/6/13 by DJ - cleaned up code
% Updated 2/26/16 by DJ - added LTO inputs/outputs

if nargin<3 || isempty(params)
    % Set parameters
    regularize = 1;
    lambda = 1e-6;
    lambdasearch = true;
    eigvalratio = 1e-4;
    vinit = zeros(size(data,1)+1,1);
    show = 0;
    LOO = false; % perform leave-one-out CV
    LTO = true; % perform leave-two-out CV
    demean = false; % subtract mean y value on each CV fold
else
    % Unpack parameters from struct
    UnpackStruct(params);
end

if ~exist('demean','var')
    demean = false;
end
if ~exist('LTO','var')
    LTO = false;
end

% Set up
[nFeats, nSamples, nTrials] = size(data);

% Set up LR on training data
truth = reshape(repmat(trialTruth(:)',nSamples,1),nTrials*nSamples,1);
trial = reshape(repmat(1:nTrials,nSamples,1),nTrials*nSamples,1);
x = data(:,:)'; % Rearrange data for logist.m [(T x trials), D]
% Perform LR
v = logist(x,truth,vinit,show,regularize,lambda,lambdasearch,eigvalratio); % weights
y = [x, ones(nSamples*nTrials,1)]*v; % classification values
% y = [x, zeros(nSamples*nTrials,1)]*v; % classification values
% Use mean y value for each trial to classify it
ymean = mean(reshape(y,nSamples,nTrials),1);
bp = bernoull(1,ymean); 
% Get Az value, ROC curve
[Az,Ry,Rx] = rocarea(bp,trialTruth); 
% fprintf('Training Az = %.2f\n',Az);
% Get fwd model
a = y \ x; % fwd model
% Create stats struct
stats.wts = v;
stats.fwdModel = a;
stats.y = y;


% Leave-One-Out Analysis (LOO)
if LOO
    for looi=1:nTrials,
        % Extract data
        xLoo = x(trial~=looi,:); % LOO data
        truthLoo = truth(trial~=looi);   % LOO truth
        % Get weights
        vLoo(:,looi) = logist(xLoo,truthLoo,vinit,show,regularize,lambda,lambdasearch,eigvalratio);
        % Get y value for left-out trial
        yLoo(:,looi) = [x(trial==looi,:), ones(nSamples,1)]*vLoo(:,looi);
        % if requested, subtract mean of training trials
        if demean
%             yAll = [xLoo, ones(size(xLoo,1),1)]*vLoo(:,looi); % Don't include left-out sample in mean
            yAll = [x, ones(size(x,1),1)]*vLoo(:,looi); % include left-out sample in mean
            yLoo(:,looi) = yLoo(:,looi) - mean(yAll);
        end

        % Use mean y value for trial to classify it
        ymean(looi) = mean(yLoo(:,looi));
        bploomean(looi)=bernoull(1,ymean(looi));
        % Get forward model
        aLoo(:,looi) = yLoo(:,looi) \ x(trial==looi,:);
    end
    % Get Az value, ROC curve
    [AzLoo,Ryloo,Rxloo] = rocarea(bploomean,trialTruth);
    % Add LOO stuff to stats struct
    stats.wtsLoo = vLoo;
    stats.fwdModelLoo = aLoo;
    stats.yLoo = yLoo;
    stats.AzLoo = AzLoo;
else
    AzLoo = NaN;    
    stats.wtsLoo = NaN;
    stats.fwdModelLoo = NaN;
    stats.yLoo = NaN;
    stats.AzLoo = NaN;
end

% Perform leave-two-out analysis
if LTO 
    % Crop to even number of trials in each class (first chronologically)    
    zeroTrials = find(truthTrials<=0);
    oneTrials = find(truthTrials<=1);
    nPerClass = min([numel(zeroTrials),numel(oneTrials)]);
    isOkTrial = ismember(trial, [zeroTrials(1:nPerClass), oneTrials(1:nPerClass)]);
    vLto = nan(nFeats,nPerClass);
    yLto = nan(nSamples,nTrials);
    bpltomean = nan(nSamples,nTrials);
    aLto = nan(nSamples,nFeats);
    for ltoi=1:nPerClass
        % Extract data
        xLto = x((trial~=zeroTrials(ltoi) & trial~=oneTrials(ltoi) & isOkTrial), :); % LTO data
        truthLto = truth((trial~=zeroTrials(ltoi) & trial~=oneTrials(ltoi) & isOkTrial), :); % LTO truth
        % Get weights
        vLto(:,ltoi) = logist(xLto,truthLto,vinit,show,regularize,lambda,lambdasearch,eigvalratio);
        % Get y value for left-out trials
        yLto(:,zeroTrials(ltoi)) = [x(trial==zeroTrials(ltoi),:), ones(nSamples,1)]*vLto(:,ltoi);
        yLto(:,oneTrials(ltoi)) = [x(trial==oneTrials(ltoi),:), ones(nSamples,1)]*vLto(:,ltoi);
        % if requested, subtract mean of training trials
        if demean
%             yAll = [xLto, ones(size(xLto,1),1)]*vLto(:,ltoi); % Don't include left-out sample in mean
            yAll = [x(isOkTrial,:), ones(size(x,1),1)]*vLto(:,ltoi); % include left-out samples in mean
            yLto(:,zeroTrials(ltoi)) = yLto(:,zeroTrials(ltoi)) - mean(yAll);
            yLto(:,oneTrials(ltoi)) = yLto(:,oneTrials(ltoi)) - mean(yAll);
        end

        % Use mean y value for trial to classify it
        ymean(zeroTrials(ltoi)) = mean(yLto(:,zeroTrials(ltoi)));
        ymean(oneTrials(ltoi)) = mean(yLto(:,oneTrials(ltoi)));        
        bpltomean(zeroTrials(ltoi))=bernoull(1,ymean(zeroTrials(ltoi)));
        bpltomean(oneTrials(ltoi))=bernoull(1,ymean(oneTrials(ltoi)));
        % Get forward model
        aLto(:,zeroTrials(ltoi)) = yLto(:,zeroTrials(ltoi)) \ x(trial==zeroTrials(ltoi),:);
        aLto(:,oneTrials(ltoi)) = yLto(:,oneTrials(ltoi)) \ x(trial==oneTrials(ltoi),:);
    end
    % Get Az value, ROC curve
    [AzLto,Rylto,Rxlto] = rocarea(bpltomean(~isnan(bpltomean)),trialTruth(~isnan(bpltomean)));
    % Add LTO stuff to stats struct
    stats.wtsLto = vLto;
    stats.fwdModelLto = aLto;
    stats.yLto = yLto;
    stats.AzLto = AzLto;
else
    AzLto = NaN;    
    stats.wtsLto = NaN;
    stats.fwdModelLto = NaN;
    stats.yLto = NaN;
    stats.AzLto = NaN;
end
