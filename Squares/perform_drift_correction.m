function [sq, drift, rotation] = perform_drift_correction(sq,outlier_cutoff,doplot)

% Compensate for eye position offset drift by re-offsetting data so that
% closest fixation lands directly on the fixation cross.
%
% [sq, drift,rotation] = perform_drift_correction(sq)
%
% INPUTS:
% - sq is a squares data struct imported using import_squares_data.
% - outlier_cutoff is a scalar indicating the number of std dev's from the
% mean that a drift/rotation value must be to be considered an outlier.
% - doplot is a binary value indicating whether you want to plot the drift
% in the current figure.
%
% OUTPUTS:
% - sq is the same data struct with fixation and saccade positions on each
% trial offset to counter the drift.
% - drift is an nx2 matrix, where n is the number of trials in sq.  sq(i,:)
% is the x and y position of the fixation cross minus the position of the
% nearest fixation to the cross on that trial.
% - rotation is an n-element vector indicating the rotation detected on
% each trial.
%
% Created 10/19/11 by DJ.
% Updated 11/16/11 by DJ - comments.
% Updated 4/20/12 by DJ - added circle anchoring point and rotation.
% Updated 4/23/12 by DJ - added legacy version for squares_v1pt4, switched  
%   to 3*stddev as outlier cutoff
% Updated 4/30/12 by DJ - added plot capability, doplot input
% Updated 5/30/12 by DJ - added outlier_cutoff input
% Updated 7/23/12 by DJ - comments.
% Updated 3/4/13 by DJ - added SquaresFix_v1pt2 support.
% Updated 2/20/14 by DJ - added SquaresFix3_v1pt3 support.

if nargin<2 || isempty(outlier_cutoff)
    outlier_cutoff = 3; % std dev's from mean required to negate a drift measure
end
if nargin<2 || isempty(doplot)
    doplot = 0;
end


Constants = GetSquaresConstants();
fix_error = nan(length(sq.trial.start_time),2);
circle_error = nan(length(sq.trial.start_time),2);
drift = zeros(length(sq.trial.start_time),2);
rotation = zeros(length(sq.trial.start_time),1);
center = [mean([Constants.LEFTCROSS_X,Constants.RIGHTCROSS_X]), ...
    mean([Constants.LEFTCROSS_Y, Constants.RIGHTCROSS_Y])]; % center of screen, used for rotation anchoring
nTrials = numel(sq.trial.start_time);

disp('Correcting drift and rotation...')
switch sq.experiment
    case {'SquaresFix_v1pt2', 'SquaresFix3_v1pt0'}
        disp('(SquaresFix Experiment: no rotation correction.)');
        for i=1:nTrials
            thisCrossTime = sq.trial.fix_time(i);
            thisStartTime= sq.trial.start_time(i);
            fix_pos = center;            
           
            % Check for fixations during fixation cross and end circle times
        %     iFix = find(sq.fixation.end_time>=thisCrossTime & sq.fixation.end_time<nextCrossTime);
            iFixStart = find(sq.fixation.end_time>=thisCrossTime & sq.fixation.start_time<thisStartTime); % all fixations that overlap with the cross period
            
            % If these fixations are not found, use the ones from the trial before
            if isempty(iFixStart)
                fprintf('Trial %d drift not found! Using previous trial''s drift...\n',i);
                if i==1 % there is no trial before, so assume no drift.
                    drift(i,:) = 0;
                else  % use the drift and rotation from the previous trial.
                    drift(i,:) = drift(i-1,:);
                end 
            else % Cross and Circle fixations are both found
                % Find closest fixations
                [~, iClosestToFix] = min(sum((repmat(fix_pos,length(iFixStart),1) - sq.fixation.position(iFixStart,:)).^2,2)); % find closest fixation to cross
                
                % Get drift and rotation
                drift(i,:) = fix_pos - sq.fixation.position(iFixStart(iClosestToFix),:);                
            end
            rotation(i,:) = 0;
        end
    
    case 'Squares_v1pt4' % OLD VERSION - DRIFT CALCULATED FROM CROSS ONLY
        disp('(Legacy version: no rotation correction.)')
        for i=1:nTrials
            thisCrossTime = sq.trial.fix_time(i);
            thisStartTime= sq.trial.start_time(i);
            if sq.trial.is_right_cross(i)
                fix_pos = [Constants.RIGHTCROSS_X, Constants.RIGHTCROSS_Y];
            else
                fix_pos = [Constants.LEFTCROSS_X, Constants.LEFTCROSS_Y];
            end

            % Check for fixations during fixation cross and end circle times
        %     iFix = find(sq.fixation.end_time>=thisCrossTime & sq.fixation.end_time<nextCrossTime);
            iFixStart = find(sq.fixation.end_time>=thisCrossTime & sq.fixation.start_time<thisStartTime); % all fixations that overlap with the cross period
            
            % If these fixations are not found, use the ones from the trial before
            if isempty(iFixStart)
                fprintf('Trial %d drift not found! Using previous trial''s drift...\n',i);
                if i==1 % there is no trial before, so assume no drift.
                    drift(i,:) = 0;
                else  % use the drift and rotation from the previous trial.
                    drift(i,:) = drift(i-1,:);
                end 
            else % Cross and Circle fixations are both found
                % Find closest fixations
                [~, iClosestToFix] = min(sum((repmat(fix_pos,length(iFixStart),1) - sq.fixation.position(iFixStart,:)).^2,2)); % find closest fixation to cross
                
                % Get drift and rotation
                drift(i,:) = fix_pos - sq.fixation.position(iFixStart(iClosestToFix),:);                
            end
            rotation(i,:) = 0;
        end
        
    otherwise % CURRENT VERSION - USE BOTH CROSS AND CIRCLE AS ANCHORS
        for i=1:nTrials
            thisCrossTime = sq.trial.fix_time(i);
            thisStartTime= sq.trial.start_time(i);
            thisEndTime = sq.trial.end_time(i);
            thisCircleTime = sq.trial.circle_time(i);
            if sq.trial.is_right_cross(i)
                fix_pos = [Constants.RIGHTCROSS_X, Constants.RIGHTCROSS_Y];
                circle_pos = [Constants.LEFTCROSS_X, Constants.LEFTCROSS_Y];
            else
                fix_pos = [Constants.LEFTCROSS_X, Constants.LEFTCROSS_Y];
                circle_pos = [Constants.RIGHTCROSS_X, Constants.RIGHTCROSS_Y];
            end

            % Check for fixations during fixation cross and end circle times
        %     iFix = find(sq.fixation.end_time>=thisCrossTime & sq.fixation.end_time<nextCrossTime);
            iFixStart = find(sq.fixation.end_time>=thisCrossTime & sq.fixation.start_time<thisStartTime); % all fixations that overlap with the cross period
%             iFixStart = find(sq.fixation.start_time<thisStartTime,1,'last'); % all fixations that overlap with the cross period
            iFixEnd = find(sq.fixation.end_time>=thisEndTime & sq.fixation.start_time<thisCircleTime); % all fixations that overlap with the circle period
%             iFixEnd = find(sq.fixation.start_time<thisCircleTime,1,'last'); % all fixations that overlap with the circle period

            % If these fixations are not found, use the ones from the trial before
            if isempty(iFixStart) || isempty(iFixEnd)
                fprintf('Trial %d drift not found! Using previous trial''s drift...\n',i);
                if i==1 % there is no trial before, so assume no drift.
                    drift(i,:) = 0;
                    rotation(i,:) = 0;
                else  % use the drift and rotation from the previous trial.
                    drift(i,:) = drift(i-1,:);
                    rotation(i,:) = rotation(i-1,:);
                end 
            else % Cross and Circle fixations are both found
                % Find closest fixations
                [~, iClosestToFix] = min(sum((repmat(fix_pos,length(iFixStart),1) - sq.fixation.position(iFixStart,:)).^2,2)); % find closest fixation to cross
                [~, iClosestToCircle] = min(sum((repmat(circle_pos,length(iFixEnd),1) - sq.fixation.position(iFixEnd,:)).^2,2)); % find closest fixation to circle

                % Get drift and rotation
                fix_error(i,:) = fix_pos - sq.fixation.position(iFixStart(iClosestToFix),:);
                circle_error(i,:) = circle_pos - sq.fixation.position(iFixEnd(iClosestToCircle),:);
                drift(i,:) = mean([fix_error(i,:); circle_error(i,:)],1);
                rotation(i,:) = atan( (sq.fixation.position(iFixStart(iClosestToFix),2)-sq.fixation.position(iFixEnd(iClosestToCircle),2)) / ...
                    (sq.fixation.position(iFixStart(iClosestToFix),1)-sq.fixation.position(iFixEnd(iClosestToCircle),1)) );
            end

        end
end

% Correct outliers
mean_drift = mean(drift,1);
mean_rotation = mean(rotation,1);
std_drift = std(drift,[],1);
std_rotation = std(rotation,[],1);
outliers = find(drift(:,1)>mean_drift(1)+outlier_cutoff*std_drift(1) | ...
    drift(:,1)<mean_drift(1)-outlier_cutoff*std_drift(1) | ...
    drift(:,2)>mean_drift(2)+outlier_cutoff*std_drift(2) | ...
    drift(:,2)<mean_drift(2)-outlier_cutoff*std_drift(2) | ...
    rotation>mean_rotation+outlier_cutoff*std_rotation | ...
    rotation<mean_rotation-outlier_cutoff*std_rotation); % >3 std away from the mean
if ~isempty(outliers)
    fprintf('Correcting %d outliers...\n',numel(outliers))
end

% Plot if requested
if doplot
    clf;
    subplot(2,1,1)
    scatter(drift(:,1),drift(:,2),'b.');
    hold on
    scatter(drift(outliers,1),drift(outliers,2),'r.');
    if ~isinf(outlier_cutoff)
        rectangle('Position',[mean_drift(1)-outlier_cutoff*std_drift(1),mean_drift(2)-outlier_cutoff*std_drift(2),...
            2*outlier_cutoff*std_drift(1),2*outlier_cutoff*std_drift(2)],'EdgeColor','r','Linestyle',':');
    end
    xlabel('x drift (pixels)')
    ylabel('y drift (pixels)')
    title(show_symbols(sprintf('%s trial drift',sq.eyeFilename)))
    
    subplot(2,1,2)
    [AX,H1,H2] = plotyy(1:nTrials, drift, 1:nTrials,rotation);
    hold(AX(1),'on'); hold(AX(2),'on');
    set(H1,'Marker','.'); set(H2,'Marker','.');
%     plot(drift(:,1),'c.-');
%     hold on
%     plot(drift(:,2),'m.-');
    xlabel('trial')
    set(get(AX(1),'Ylabel'),'String','Drift (pixels)')
    set(get(AX(2),'Ylabel'),'String','Rotation (radians)')
end

% Correct outliers
nonoutliers = setdiff(1:nTrials,outliers);
drift(outliers,:) = interp1(nonoutliers,drift(nonoutliers,:),outliers,'linear','extrap');
rotation(outliers,:) = interp1(nonoutliers,rotation(nonoutliers,:),outliers,'linear','extrap');

% for i=1:numel(outliers)
%     fprintf('Trial %d...\n',outliers(i))    
%     if outliers(i)>1 % use the drift and rotation of the previous trial, if possible
%         drift(outliers(i),:) = drift(outliers(i)-1,:);
%         rotation(outliers(i),:) = rotation(outliers(i)-1,:);   
%     elseif all(outliers~=outliers(i)+1) % use the drift and rotation of the following trial, if it's not an outlier
%         drift(outliers(i),:) = drift(outliers(i)+1,:);
%         rotation(outliers(i),:) = rotation(outliers(i)+1,:);      
%     else % set the drift and rotation to 0 (as a last resort).
%         drift(outliers(i),:) = 0;
%         rotation(outliers(i),:) = 0;
%     end
%         
% end

% Plot corrected drifts
if doplot
%     subplot(2,1,2)
    plot(AX(1),drift,'o');
    plot(AX(2),rotation,'o');   
%     plot(AX(1),drift(:,1),'co');
%     plot(AX(1),drift(:,2),'mo');
    if ~isinf(outlier_cutoff)
        plot(AX(1),[0 length(drift)+1],[mean_drift(1)-outlier_cutoff*std_drift(1),mean_drift(1)-outlier_cutoff*std_drift(1)],'b--');
        plot(AX(1),[0 length(drift)+1],[mean_drift(1)+outlier_cutoff*std_drift(1),mean_drift(1)+outlier_cutoff*std_drift(1)],'b--');
        plot(AX(1),[0 length(drift)+1],[mean_drift(2)-outlier_cutoff*std_drift(2),mean_drift(2)-outlier_cutoff*std_drift(2)],'g--');
        plot(AX(1),[0 length(drift)+1],[mean_drift(2)+outlier_cutoff*std_drift(2),mean_drift(2)+outlier_cutoff*std_drift(2)],'g--');
        plot(AX(2),[0 length(rotation)+1],[mean_rotation-outlier_cutoff*std_rotation,mean_rotation-outlier_cutoff*std_rotation],'r--');
        plot(AX(2),[0 length(rotation)+1],[mean_rotation+outlier_cutoff*std_rotation,mean_rotation+outlier_cutoff*std_rotation],'r--');
    end
    legend('x','y','rotation')  
end

% Apply drift correction to saccades and fixations
for i=1:nTrials
    thisCrossTime = sq.trial.fix_time(i);
    if i==numel(sq.trial.start_time)
        nextCrossTime = Inf;
    else
        nextCrossTime = sq.trial.fix_time(i+1);
    end
    % Use drift to update all postions
    iSac = find(sq.saccade.start_time>=thisCrossTime & sq.saccade.start_time<nextCrossTime);
    sq.saccade.start_position(iSac,:) = sq.saccade.start_position(iSac,:) + repmat(drift(i,:),numel(iSac),1);
    sq.saccade.end_position(iSac,:) = sq.saccade.end_position(iSac,:) + repmat(drift(i,:),numel(iSac),1);
%     iFix = find(sq.fixation_start_time>=thisCrossTime & sq.fixation_start_time<nextCrossTime);
    sq.fixation.position(iSac,:) = sq.fixation.position(iSac,:) + repmat(drift(i,:),numel(iSac),1);
    
    % Use rotation to update all positions about the center of the screen
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
end





disp('Success!')