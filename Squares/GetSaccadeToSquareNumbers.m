function iSaccades = GetSaccadeToSquareNumbers(y)

% iSaccades = GetSaccadeToSquareNumbers(y)
%
% INPUTS:
% -y is an n-element vector of structs, where y(i) is the behavioral data 
% for session i.
%
% OUTPUTS:
% -iSaccades is an n-element cell vector, and iSaccades{i} is an
% nTrials x nSquares matrix of saccade indices in y(i).saccade.  If 
% iSaccades{i}(j,k) = l, then saccade #l is the first saccade in session i 
% for which trialnum==j and squarenum==k.
%
% Created 6/18/12 by DJ to help find matching times for LogOddsRatio events

% Set up
nSessions = length(y);
iSaccades = cell(1,nSessions);

% Find saccades
for i=1:nSessions
    [nTrials, nSquares] = size(y(i).trial.is_target_color); % get dimensions
    iSaccades{i} = nan(nTrials,nSquares); % any saccade not found will be a nan
    for j=1:nTrials
        for k=1:nSquares
            iSac = find(y(i).saccade.trialnum==j & y(i).saccade.squarenum==k,1); % get just the first one to each square
            if ~isempty(iSac)
                iSaccades{i}(j,k) = iSac; % add to matrix
            end
        end
    end
end