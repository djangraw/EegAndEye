function sq = undo_drift_correction(sq)

% Reverse the effects of perform_drift_correction.
%
% sq = undo_drift_correction(sq)
%
% INPUTS:
% - sq is a squares dataset whose drift correction you want to undo. The
% drift parameters should be found in x.trial.drift and x.trial.rotation.
% 
% OUTPUTS:
% - sq is the input dataset with drift correction undone. x.trial.drift and
% x.trial.rotation are set to zeros.
%
% Created 2/8/13 by DJ.

% Get constants
nTrials = length(sq.trial.start_time);
Constants = GetSquaresConstants();
center = [mean([Constants.LEFTCROSS_X,Constants.RIGHTCROSS_X]), ...
    mean([Constants.LEFTCROSS_Y, Constants.RIGHTCROSS_Y])]; % center of screen, used for rotation anchoring

% Parse out drift & rotation
drift_toundo = sq.trial.drift;
if isfield(sq.trial,'rotation')
    rotation_toundo = sq.trial.rotation;
else
    rotation_toundo = zeros(nTrials,1);
end

% reverse drift & rotation
drift = -drift_toundo;
rotation = -rotation_toundo;

% Apply drift correction to saccades and fixations
for i=1:nTrials
    % Find saccades from this trial
    thisCrossTime = sq.trial.fix_time(i);
    if i==numel(sq.trial.start_time)
        nextCrossTime = Inf;
    else
        nextCrossTime = sq.trial.fix_time(i+1);
    end
    iSac = find(sq.saccade.start_time>=thisCrossTime & sq.saccade.start_time<nextCrossTime);

    % perform_drift_correction applies (1)drift and (2)rotation, so to undo
    % it we reverse the operations.
    
    % (1) Use rotation to update all positions about the center of the screen
    rotMatrix = [cos(rotation(i)), -sin(rotation(i)); sin(rotation(i)), cos(rotation(i))];
    
    sq.saccade.start_position(iSac,:) = sq.saccade.start_position(iSac,:) - repmat(center,numel(iSac),1);
    sq.saccade.end_position(iSac,:) = sq.saccade.end_position(iSac,:) - repmat(center,numel(iSac),1);
    sq.fixation.position(iSac,:) = sq.fixation.position(iSac,:) - repmat(center,numel(iSac),1);
    
    sq.saccade.start_position(iSac,:) = sq.saccade.start_position(iSac,:) * rotMatrix;
    sq.saccade.end_position(iSac,:) = sq.saccade.end_position(iSac,:) * rotMatrix;
    sq.fixation.position(iSac,:) = sq.fixation.position(iSac,:) * rotMatrix;
    
    sq.saccade.start_position(iSac,:) = sq.saccade.start_position(iSac,:) + repmat(center,numel(iSac),1);
    sq.saccade.end_position(iSac,:) = sq.saccade.end_position(iSac,:) + repmat(center,numel(iSac),1);
    sq.fixation.position(iSac,:) = sq.fixation.position(iSac,:) + repmat(center,numel(iSac),1);
    
    % (2) Use drift to update all postions
    sq.saccade.start_position(iSac,:) = sq.saccade.start_position(iSac,:) + repmat(drift(i,:),numel(iSac),1);
    sq.saccade.end_position(iSac,:) = sq.saccade.end_position(iSac,:) + repmat(drift(i,:),numel(iSac),1);
%     iFix = find(sq.fixation_start_time>=thisCrossTime & sq.fixation_start_time<nextCrossTime);
    sq.fixation.position(iSac,:) = sq.fixation.position(iSac,:) + repmat(drift(i,:),numel(iSac),1);

end

% zero out data struct's drift/rotation fields
sq.trial.drift = zeros(size(drift));
sq.trial.rotation = zeros(size(rotation));