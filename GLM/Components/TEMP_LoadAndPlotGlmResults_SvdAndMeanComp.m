% TEMP_LoadAndPlotGlmResults_SvdAndMeanComp.m
%
% Loads GLM results from all subjects, plots the group SVD results, and
% plots the group results for a "mean component" where the weight is
% constant across all electrodes, a sort of grand average.
%
% Created 2/22/13 by DJ for one-time use.
% Updated 4/15/13 by DJ - added prefix for sf use
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 3/19/14 by DJ - tResponse in cells


%%% Squares %%%
% prefix = 'sq';
% subjects = [9:11 13:15 17:19];
% glmType = 'LREvents-SqNum-NewType-v2pt3';

%%% SquaresFix %%%
prefix = 'sf';
subjects = [5];
glmType = 'SqNum-Type-v2pt1';

%%
clear R
for i=1:numel(subjects)
    R(i) = load(sprintf('%s-%d-GLMresults-%s.mat',prefix,subjects(i),glmType));
end

%% Plot SVD results
[wts, tc,~,~,~,tResponse] = GetGroupSvdResults(R,3,size(R(1).responseFns{iLevel},3),[0 500]);
MakeFigureTitle(sprintf('%s, %d subjects, [0 500] ms SVD',glmType,numel(subjects)));

%% stats for svd
GetCrossSubjectStats(tc,tResponse,R(1).regressor_events{R(1).iLevel});
MakeFigureTitle(sprintf('%s, %d subjects, [0 500] ms SVD',glmType,numel(subjects)));

%% Get average across all electrodes

icaweights = ones(1,81)/81;
icasphere = eye(81);
icawinv = pinv(icaweights*icasphere);

scaling = repmat(sqrt(sum(icawinv(:,:).^2))', [1 size(icaweights,2)]);
    icaweights = icaweights .* scaling;
    icawinv = pinv(icaweights * icasphere);


clear I
for i=1:numel(subjects)
    I(i).icaweights = icaweights;
    I(i).icasphere = icasphere;
    I(i).icawinv = icawinv;
end

[wts, tc] = GetGroupIcaResults(R,1,6,I);

MakeFigureTitle(sprintf('%s, %d subjects, mean component',glmType,numel(subjects)));

%% stats for avg component
GetCrossSubjectStats(tc,tResponse,R(1).regressor_events{R(1).iLevel});
% MakeFigureTitle(sprintf('%s, %d subjects, mean component',glmType,numel(subjects)));

%% Plot ERP Grid
RF = [];
for i=1:numel(R)
    RF = cat(4,R,R(i).responseFns{R(i).iLevel});
end
PlotResponseFnsGrid(RF,R(1).regressor_events{R(1).iLevel},tResponse,R(1).EEG.chanlocs,{'F3' 'FZ' 'F4'; 'C3' 'CZ' 'C4'; 'P3' 'PZ' 'P4'; 'O1' 'OZ' 'O2'});
