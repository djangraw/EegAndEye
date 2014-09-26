% ProduceIeeeSaccadeTimes
%
% Plots the histograms of times at which the subject made any saccade in 
% target and distractor trials, relative to the t=0 of the epoch.
% This plot was Figure 5 in DJ's 2011 IEEE/EMBS NE conference submission.
%
% Created 3/15/11 by DJ 

%% Set up
subjects = [6 2 7];

% Load data
LoadAllEpochs; % alter this program to get stimulus-locked or saccade-locked analysis.

% Make plot
Figure;
GetEyeAverages(subjects,[0 1 0 1 0 1]);
