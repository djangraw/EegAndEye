function contrastMatrix = MakeContrastMatrix_flex(events,event_list,event_weights,tResponse,baseline_win)

% Creates a matrix to multiply by 'beta' values to combine various events.
%
% contrastMatrix =
% MakeContrastMatrix_flex(events,event_list,event_weights,tResponse,baseline_win)
%
% INPUTS:
% - events is an N-element cell array of strings indicating the event type
% of each response function in the array this contrast matrix is built for.
% - event_list is an M-element cell array of strings indicating the
% event(s) that should be added together in this contrast.
% - event_weights is a P x M matrix in which each row event_weights(i,:) 
% indicates the weight that should be assigned to the corresponding events
% in event_list in contrast i.
% - tResponse is a T-element vector indicating the times of each sample in
% the response functions.
% - baseline_win is a 2-element vector indicating the time window that
% should be considered baseline and subtracted from each event
% participating in the contrast.
%
% OUTPUTS:
% - contrastMatrix is a matrix of size [(T*N),(T*P)] that can be used by
% programs like GetGlmZscore.
%
% Created 4/9/13 by DJ.
% Updated 4/22/13 by DJ - 1:T, not 0:T-1; and jBaseline, not jBaseline-1.
% Updated 3/21/14 by DJ - _flex version

% Handle inputs and defaults
if ischar(events)
    events = {events}; % convert from string to cell array
end
if ischar(event_list)
    event_list = {event_list}; % convert from string to cell array
end
if nargin<3 || isempty(event_weights)
    event_weights = ones(size(event_list));
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
P = size(event_weights,1); % # contrasts

% get indices of events and baselines
kEvents = nan(1,numel(event_list));
for i=1:numel(event_list)
    kEvents(i) = find(strcmp(event_list{i},events));
end
% kEvents = find(ismember(events,event_list));
if numel(kEvents) ~= numel(event_list) % Make sure all events were found
    error('some events not found!')
end   
jBaseline = find(tResponse>=baseline_win(1) & tResponse<=baseline_win(2));

% Get baseline weights
baseline_weight = 1/numel(jBaseline);


% Set up contrast matrix
contrastMatrix = zeros(T*N,T*P);
for i = 1:P % for each contrast
    % Determine which columns to fill
    iCols = (i-1)*T + (1:T);
    % Insert blocks of contrast matrix
    for j=1:numel(kEvents)
        contrastMatrix((kEvents(j)-1)*T + (1:T),iCols) = event_weights(i,j)*eye(T);
        contrastMatrix((kEvents(j)-1)*T + jBaseline,iCols) = contrastMatrix((kEvents(j)-1)*T + jBaseline,:) - event_weights(i,j)*baseline_weight;
    end
end
