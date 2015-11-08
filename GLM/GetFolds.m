function [iTrain,iTest] = GetFolds(n,nFolds,useRandomFolds,nRandomizations)

% Get the training and testing trials for each cross-validation fold.
% [iTrain,iTest] = GetFolds(n,nFolds,useRandomFolds,nRandomizations)
%
% INPUTS:
% -n is a scalar indicating the number of trials in your dataset.
% -nFolds is a scalar indicating the desired number of cross-validation 
% folds.
% -useRandomFolds is a binary value indicating whether the trials in each
% fold should be selected randomly (1) or grouped chronologically (0).
% -nRandomizations is a scalar indicating the number of times the data
% should be split into folds.
%
% OUTPUTS:
% -iTrain is an [n x nRandomizations] array of cells in which iTrain{i,j}
% contains a list of the training trials for fold i in randomization j.
% -iTest is an [n x nRandomizations] array of cells in which iTest{i,j}
% contains a list of the testing trials for fold i in randomization j.
%
% Created 2/27/15 by DJ.

% Declare defaults
if ~exist('useRandomFolds','var') || isempty(useRandomFolds)
    useRandomFolds = false;
end
if ~exist('nRandomizations','var') || isempty(nRandomizations)
    nRandomizations = 1;
end

% set up
iTest = cell(nFolds,nRandomizations);
iTrain = cell(nFolds,nRandomizations);
folds = ceil(linspace(eps,nFolds,n)); % which fold should each trial be in?

for j=1:nRandomizations
    % if randomizing, scramble the folds
    if useRandomFolds    
        thisFolds = folds(randperm(numel(folds)));
    else % otherwise, use them as is
        thisFolds = folds;
    end
    % get list of train/test indices
    for i=1:nFolds
        iTest{i,j} = find(thisFolds==i);        
        iTrain{i,j} = setdiff(1:n,iTest{i,j});
    end
end