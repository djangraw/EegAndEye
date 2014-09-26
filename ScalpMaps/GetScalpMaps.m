function [results, meanResults] = GetScalpMaps(responses,tResponses,tPlots,binWidth,cthresh)

% Gets the mean responses within a set of time bins (for easy input into
% PlotScalpMaps.m).
%
% [results, meanResults] =
% GetScalpMaps(responses,tResponses,tPlots,binWidth,cthresh)
%
% INPUTS:
% - responses is a DxTxN matrix of data, where D = # channels T = # time 
%   points, N = # conditions.
% - tResponses is a T-element vector indicating the time of each sample. It
%   will be used to index the time bins [default: 1:size(responses,2)].
% - tPlots is an M-element vector indicating the center of the bins you're
%   interested in (in the same units as tResponses). 
%   [default: mean(tResponses)]
% - binWidth is a scalar indicating the width (in the same units as
%   tResponses). [default: inf (average all time points)]
% - cthresh is a scalar indicating the cutoff threshold you'd like to use.
%   All results whose magnitudes are below this threshold will be set to 0.
%   [default: 0 (no cutoff)]
%
% OUTPUTS:
% - results is a DxMxN matrix of mean responses in each time window.
% - meanResults is the average across conditions. This average is taken
%   before the cutoff threshold is applied, and the cutoff threshold is
%   then applied to meanResults as well.
%
% Created 4/29/13 by DJ.
% Updated 5/23/13 by DJ - comments.

% Handle defaults
if nargin<2 || isempty(tResponses) % if time of samples is unspecified
    tResponses = 1:size(responses,2); % use sample indices
end
if nargin<3 || isempty(tPlots) % if time of plots is unspecified
    tPlots = mean(tResponses); % use mean time of results
end
if nargin<4 || isempty(binWidth) % if bin width is not specified
    if numel(tPlots)==1
        binWidth = inf; % average all data points
    else
        binWidth = median(diff(tPlots)); % use all data between plot times
    end
end
if nargin<5 || isempty(cthresh)
    cthresh = 0; % no cutoff threshold
end

% Set up
[D,T,N] = size(responses); % # elecs, # time points, # subjects/conditions
M = numel(tPlots); % # plots
results = nan(D,M,N);

% Take average in each window
for i=1:M % window
    tWin = [tPlots(i)-binWidth/2, tPlots(i)+binWidth/2];
    isInWin = tResponses>=tWin(1) & tResponses<tWin(2);
    results(:,i,:) = mean(responses(:,isInWin,:),2);        
end

% Take mean across conditions
meanResults = mean(results,3);

% threshold results
results(abs(results)<cthresh) = 0;
meanResults(abs(meanResults)<cthresh) = 0;