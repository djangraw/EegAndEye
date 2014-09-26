% Test LOR change for each stimulus type
%
% Created 3/26/14 by DJ.

% Declare stim types of interest
types = {'D0','T0','D1','T1','D2'};
isTarget = [0 1 0 1 0];
nTargetsSoFar = [0 1 1 2 2];

% Make fake dataset with all trial types
nSquares = 5;
nTargetsForComplete = 2;
% Get all the possible sets of target colors
colors = [1;2;3];
for i=2:nSquares
    newCol = [ones(size(colors,1),1); repmat(2,size(colors,1),1); repmat(3,size(colors,1),1)];
    colors = cat(2,newCol,repmat(colors,3,1));
end
y = [];
y.trial.is_target_color = (colors==1);
y.trial.is_target_trial = sum(y.trial.is_target_color,2)>=nTargetsForComplete;
y.trial.target_squares_sofar = cumsum(y.trial.is_target_color,2);
y.trial.is_right_cross = false(size(colors,1),1);


% Get log odds ratio and other trial info
[logOddsBefore,logOddsAfter,logOddsChange] = GetLogOddsRatio(y);
allLOC = cat(1,logOddsChange{:});
allLOA = cat(1,logOddsAfter{:});
allIsTarg = []; % is target
allTSF = []; % targets so far
for i=1:numel(y)
    allIsTarg = [allIsTarg; y(i).trial.is_target_color];
    allTSF = [allTSF; y(i).trial.target_squares_sofar];
end
allIsDecided = isinf(allLOA); % can decision be made after this stimulus

% Find the mean and distribution of log odds changes for each stimulus type
nTrials = length(y(1).trial.is_target_trial);
for j=1:numel(types)
    isThisType = (allIsTarg==isTarget(j) & allTSF==nTargetsSoFar(j) & ~allIsDecided);
    fprintf('%s: %d\n',types{j},sum(isThisType(:)));
    LOCvalues = unique(allLOC(isThisType));
    fprintf('---type %s: mean LOC = %.3f\n',types{j},mean(allLOC(isThisType)));
    for k=1:numel(LOCvalues)
        fprintf('LOC =  %.3f: %.1f%% trials\n',LOCvalues(k),numel(allLOC(isThisType & allLOC==LOCvalues(k)))/nTrials*100);
    end
end