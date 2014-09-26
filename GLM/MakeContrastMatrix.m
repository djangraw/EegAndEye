function contrastMatrix = MakeContrastMatrix(events,events_plus,events_minus,tResponse,baseline_win)

% Creates a matrix to multiply by 'beta' values to combine various events.
%
% contrastMatrix =
% MakeContrastMatrix(events,events_plus,events_minus,tResponse,baseline_win)
%
% INPUTS:
% - events is an N-element cell array of strings indicating the event type
% of each 
% - events_plus is a string or cell array of strings indicating the
% event(s) that should be added together in this contrast.
% - events_minus is a string or cell array of strings indicating the
% event(s) that should be subtracted from this contrast.
% - tResponse is a T-element vector indicating the times of each sample in
% the response functions.
% - baseline_win is a 2-element vector indicating the time window that
% should be considered baseline and subtracted from each event
% participating in the contrast.
%
% OUTPUTS:
% - contrastMatrix is a matrix of size [(T*N),T] that can be used by
% programs like GetGlmZscore.
%
% Created 4/9/13 by DJ.
% Updated 4/22/13 by DJ - 1:T, not 0:T-1; and jBaseline, not jBaseline-1.

% Handle inputs and defaults
if ischar(events)
    events = {events}; % convert from string to cell array
end
if ischar(events_plus)
    events_plus = {events_plus}; % convert from string to cell array
end
if ischar(events_minus)
    events_minus = {events_minus}; % convert from string to cell array
end
if nargin<4 || isempty(tResponse)
    tResponse = 0; % one-element contrast
end
if nargin<5 || isempty(baseline_win)
    baseline_win = [0 -1]; % no baseline subtraction
end

% Infer size of contrast matrix
T = numel(tResponse); % # timepoints
N = numel(events); % # events

% get indices of events and baselines
kPlus = find(ismember(events,events_plus));
kMinus = find(ismember(events,events_minus));
if numel(kPlus) ~= sum(~cellfun('isempty',events_plus)) || numel(kMinus) ~= sum(~cellfun('isempty',events_minus)) % Make sure all events were found
    warning('some events not found!')
end   
jBaseline = find(tResponse>=baseline_win(1) & tResponse<=baseline_win(2));

% Get baseline weights
baseline_weight = 1/numel(jBaseline);

% Set up contrast matrix
contrastMatrix = zeros(T*N,T);
% make plus blocks
for i=1:numel(kPlus)
    contrastMatrix((kPlus(i)-1)*T + (1:T),:) = eye(T);
    contrastMatrix((kPlus(i)-1)*T + jBaseline,:) = contrastMatrix((kPlus(i)-1)*T + jBaseline,:) -baseline_weight;
end
% Make minus blocks
for i=1:numel(kMinus)
    contrastMatrix((kMinus(i)-1)*T + (1:T),:) = -eye(T);
    contrastMatrix((kMinus(i)-1)*T + jBaseline,:) = contrastMatrix((kMinus(i)-1)*T + jBaseline,:) + baseline_weight;
end