function [logOddsBefore,logOddsAfter,logOddsChange] = GetLogOddsRatio(y)

% Get log odds ratio (of being a complete trial) for every saccade in a 
% squares behavior file.
%
% logOdds = GetLogOddsRatio(y)
%
% INPUTS:
% - y is a 1xn struct array of squares behavior files, with 5 squares and 3
% possible colors in each trial.
%
% OUTPUTS:
% - logOdds is a 1xn cell array, where each cell contains an (nTrials x 5) 
% matrix where logOdds{i}(j,k) holds the log odds ratio of being a complete
% trial just before seeing square k on trial j in session i.
% 
% Created 6/15/12 by DJ.

%% GET SEQUENCES AND LOG-ODDS!
nSquares = size(y(1).trial.is_target_color,2);
% Get all the possible sets of target colors
colors = [1;2;3];
for i=2:nSquares
    newCol = [ones(size(colors,1),1); repmat(2,size(colors,1),1); repmat(3,size(colors,1),1)];
    colors = cat(2,newCol,repmat(colors,3,1));
end

% find whether each stim is a target and whether each trial is complete
nTargetsForComplete = max(sum(y(1).trial.is_target_color(~y(1).trial.is_target_trial,:),2))+1;
isTarget = colors==1;
isComplete = sum(isTarget,2)>=nTargetsForComplete;

% Find all unique sequences (of length 0:nSquares-1) and the probability that the
% trial will be complete given that sequence.
sequences = cell(1,nSquares);
seqPComplete = cell(1,nSquares);
seqLogOdds = cell(1,nSquares);

sequences{1} = [];
seqPComplete{1} = mean(isComplete);

for i=2:nSquares
    sequences{i} = unique(isTarget(:,1:i-1),'rows');
    seqPComplete{i} = nan(size(sequences{i},1),1);
    for j=1:size(sequences{i},1)
        matches = ismember(isTarget(:,1:i-1),sequences{i}(j,:),'rows');
        seqPComplete{i}(j) = mean(isComplete(matches));
    end
end
for i=1:nSquares
    seqLogOdds{i} = log(seqPComplete{i}./(1-seqPComplete{i}));
    % Display results
%     disp(num2str([sequences{i}, seqPComplete{i}, seqLogOdds{i}]))
end

%% GET TRUE SEQUENCES

nSessions = length(y);
[logOddsBefore, logOddsAfter, logOddsChange] = deal(cell(1,nSessions));

for i=1:nSessions
    truth = y(i).trial.is_target_color;
    if isfield(y(i).trial,'is_right_cross');
        is_reversed = y(i).trial.is_right_cross;
    else
        is_reversed = false(size(y(i).trial.start_time));
    end
    logOddsBefore{i}(~is_reversed,1) = seqLogOdds{1};
    logOddsBefore{i}(is_reversed,5) = seqLogOdds{1};
    for j=1:size(truth,1)             
        for k=2:size(truth,2)
            if is_reversed(j)
                logOddsBefore{i}(j,5-k+1) = seqLogOdds{k}(ismember(sequences{k},truth(j,5:-1:(5-k+2)),'rows'));
            else
                logOddsBefore{i}(j,k) = seqLogOdds{k}(ismember(sequences{k},truth(j,1:k-1),'rows'));
            end
        end
    end
    logOddsFinal = repmat(-inf,length(y(i).trial.is_target_trial),1);
    logOddsFinal(y(i).trial.is_target_trial) = inf;
    logOddsAfter{i} = [logOddsBefore{i}(:,2:end), logOddsFinal];
    logOddsChange{i} = logOddsAfter{i} - logOddsBefore{i};
end
  
    



