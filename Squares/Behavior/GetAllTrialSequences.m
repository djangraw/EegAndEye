function [trials, chance] = GetAllTrialSequences(event_types,nSquares)

% Returns a list of the possible trial sequences in the squares experiment.
%
% [trials, chance] = GetAllTrialSequences
%
% OUTPUTS:
% - trials - each row is a unique trial sequence, in which:
%   0=D-0T, 1=Integ, 2=D-1T, 3=Compl, 4=Irrel.
% - chance is a vector in which chance(i) is the probability of the 
%   sequence trials(:,i) happening.
%
% Created 3/19/12 by DJ.
% Updated 5/15/12 by DJ - made circle event types optional
% Updated 2/27/13 by DJ - added Extra event

% Handle defaults
if nargin<1
    event_types = {'Dist-0T' 'Dist-1T' 'Integ' 'Compl' 'Irrel' ...
        'Circle-D' 'Circle-T'};
end
if nargin<2
    nSquares = 5;
end

% Get all the possible sets of target colors
colors = [1;2;3];
for i=2:nSquares
    newCol = [ones(size(colors,1),1); repmat(2,size(colors,1),1); repmat(3,size(colors,1),1)];
    colors = cat(2,newCol,repmat(colors,3,1));
end

% find whether each stim is a target and how many have already been seen
isTarget = colors==1;
nTargets = cumsum(isTarget, 2);

% Turn this into "stimulus type" codes
alltrials = zeros(size(colors));
alltrials(~isTarget & nTargets==0) = find(strcmp('Dist-0T',event_types)); % D-0T
alltrials(isTarget & nTargets==1) = find(strcmp('Integ',event_types)); % Integration
alltrials(~isTarget & nTargets==1) = find(strcmp('Dist-1T',event_types)); % D-1T
alltrials(isTarget & nTargets==2) = find(strcmp('Compl',event_types)); % Completion
alltrials((~isTarget & nTargets==2) | nTargets>2) = find(strcmp('Irrel',event_types)); % Irrelevant
alltrials(isTarget & nTargets==3) = find(strcmp('Extra',event_types)); % Completion

% Get unique trial sequences
sq_trials = unique(alltrials,'rows');
% Get circle events
sq_isComplete = sum(sq_trials==find(strcmp('Compl',event_types)),2)>0;
if any(strncmp('Circle',event_types,5))
    sq_circle = repmat(find(strcmp('Circle-D',event_types)),length(sq_isComplete),1);
    sq_circle(sq_isComplete) = find(strcmp('Circle-T',event_types));
else
    sq_circle = [];
end
% Put them together
trials = [sq_trials, sq_circle];

% Find how prevalent each sequence is
nExamples = zeros(size(trials,1),1);
for i=1:size(trials,1);
    nExamples(i) = sum(ismember(alltrials,trials(i,1:nSquares), 'rows'));
end
    
chance = nExamples/sum(nExamples);