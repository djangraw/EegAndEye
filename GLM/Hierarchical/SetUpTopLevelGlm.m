function [contrastFns, contrastVar, contrastZ]  = SetUpTopLevelGlm(results,events_plus,events_minus,baseline_win,iLevel)

% [contrastFns, contrastVar, contrastZ]  =
% SetUpTopLevelGlm(results,events_plus,events_minus,baseline_win,iLevel)
%
% INPUTS:
% -results is an N-element vector of GLM results structs.
% -events_plus and events_minus are strings or M-element cell arrays of 
% strings. They indicate the events you want to add and subtract in the
% contrast of interest.
% -baseline_win is a 2-element vector indicating the start and end time (in
% ms) of the window you want to subtract from each contrast event as
% baseline.
% -iLevel is a scalar indicating what level of analysis you want to use
% (default = results(1).iLevel)
%
% OUTPUTS:
% -contrastFns is a DxTxMxN matrix inndicating the contrast function for
% each subject.
% -contrastVar is a DxTxMxN matrix inndicating the variance of each 
% contrast function estimate.
% -contrastZ is the single-subject Z score for each contrast function 
% estimate. No multiple-comparisons correction is performed.
%
% Created 4/18/13 by DJ.
% Updated 4/24/13 by DJ - contrastZ output, comments.
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 3/19/14 by DJ - tResponse in cells

if nargin<5 || isempty(iLevel)
    iLevel = results(1).iLevel;
end

results = UpdateGlmResultsFormat(results);
% Set up
N = numel(results);
D = size(results(1).responseFns{iLevel},1);
T = size(results(1).responseFns{iLevel},2);
tResponse = results(1).tResponse{iLevel};
multcompare = 'none';

% Get contrasts
contrasts = MakeContrastMatrix(results(1).regressor_events{iLevel}, events_plus, events_minus, tResponse, baseline_win);
M = size(contrasts,2)/T; % number of contrasts

% find contrast functions
contrastFns = nan(D,T,M,N);
contrastVar = nan(D,T,M,N);
contrastZ = nan(D,T,M,N);
for i=1:N
    % prepare figure
    figure(100+i); clf;
    MakeFigureTitle(sprintf('%s, %s - %s contrast',results(i).filenames{iLevel},events_plus,events_minus),1);
    % Calculate and plot
    [contrastZ(:,:,:,i),~,~,contrastFns(:,:,:,i),contrastVar(:,:,:,i)] = GetGlmZscore(results(i).EEG,results(i),multcompare,contrasts);
end
