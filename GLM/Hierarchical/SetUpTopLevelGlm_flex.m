function [contrastFns, contrastVar, contrastZ]  = SetUpTopLevelGlm_flex(results,event_list,event_weights,baseline_win,iLevel,doPlot)

% [contrastFns, contrastVar, contrastZ]  =
% SetUpTopLevelGlm_flex(results,events_plus,events_minus,baseline_win,iLevel)
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
% - doPlot is a binary value indicating whether you'd like to plot the
% results in the current figure. [default is true]
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
% Updated 3/21/14 by DJ - flex version
% Updated 8/7/14 by DJ - added doPlot input.

if nargin<4 || isempty(baseline_win)
    baseline_win = []; % use default setting
end
if nargin<5 || isempty(iLevel)
    iLevel = results(1).iLevel;
end
if nargin<6 || isempty(doPlot)
    doPlot = true;
end

results = UpdateGlmResultsFormat(results);
% Set up
N = numel(results);
D = size(results(1).responseFns{iLevel},1);
T = size(results(1).responseFns{iLevel},2);
tResponse = results(1).tResponse{iLevel};
multcompare = 'none';

% Get contrasts
contrasts = MakeContrastMatrix_flex(results(1).regressor_events{iLevel}, event_list, event_weights, tResponse, baseline_win);
M = size(contrasts,2)/T; % number of contrasts

% find contrast functions
contrastFns = nan(D,T,M,N);
contrastVar = nan(D,T,M,N);
contrastZ = nan(D,T,M,N);
for i=1:N
    fprintf('---subject %d/%d...\n',i,N);
    % prepare figure    
    if doPlot
        figure(100+i); clf;
        event_list_str = sprintf('%s, ',event_list{:});
        event_weight_str = sprintf('%.1g, ',event_weights);    
        MakeFigureTitle(sprintf('%s, [%s] = [%s] contrast',results(i).filenames{iLevel},event_list_str(1:end-2),event_weight_str(1:end-2)),1);
    end
    % Calculate and plot
    results(i).iLevel = iLevel;
    [contrastZ(:,:,:,i),~,~,contrastFns(:,:,:,i),contrastVar(:,:,:,i)] = GetGlmZscore(results(i).EEG,results(i),multcompare,contrasts,doPlot);
end
disp('Done!');