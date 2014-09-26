function [dwell, isToTarget] = GetDwellTimes(subject,sessions)

% Get the time the eye dwells on each object in a 3DS task.
%
% [dwell, isToTarget] = GetDwellTimes(subject,sessions)
% [dwell, isToTarget] = GetDwellTimes(y)
%
% INPUTS:
% -subject is a scalar indicating the subject number.
% -sessions is an n-element vector indicating the session numbers for that
% subject. The program loads files called '3DS-<subject>-<session>.mat'.
% -y is an n-element vector of 3DS data structs for one subject.
%
% OUTPUTS:
% -dwell is an n-element cell vector, where each cell contains a vector of
% the time the subject spent looking at each object in the session (in ms).
% -isToTarget is an n-element cell vector, where each cell contains a 
% boolean vector indicating whether the object being viewed was a target.
%
% Created 5/21/13 by DJ.
% Updated 6/5/13 by DJ - comments.
% Updated 8/13/13 by DJ - made DIST_CUTOFF a variable, not hard-coded

% Set options
DIST_CUTOFF = 100;
option_plot = false;
% Handle inputs
if isstruct(subject)
    y = subject;
else
    y = loadBehaviorData(subject,sessions,'3DS');
end

% Set up
if option_plot
    cla; hold on;
end
% Main loop
dwell = cell(1,numel(y));
isToTarget = cell(1,numel(y));
for i=1:numel(y)    
    % Get times when the subj started and stopped looking at object
    fprintf('Session %d/%d...\n',i,numel(y));
%     foo = load(sprintf('3DS-%d-%d-eyepos.mat',subject,sessions(i)));
%     times = (1:length(foo.eyepos))+y(i).eyelink.record_time-1;
%     dist = DistToObject(y(i),foo.eyepos,times);
    times = y(i).eyelink.saccade_times;
    dist = DistToObject(y(i),y(i).eyelink.saccade_positions,times);
    isOnObject = dist<DIST_CUTOFF;
    onTimes = times(diff([0; isOnObject])>0);
    offTimes = times(diff([0; isOnObject])<0);
    
    % Get object lifetimes
    objectOnTimes = y(i).eyelink.object_events(y(i).eyelink.object_events(:,2)<1000,1);
    objectOffTimes = y(i).eyelink.object_events(y(i).eyelink.object_events(:,2)>1000,1);    
    objectNumbers = y(i).eyelink.object_events(y(i).eyelink.object_events(:,2)>1000,2)-1000;
    objectIsTarget = strcmp('TargetObject',{y(i).objects(objectNumbers).tag});
    for j=1:length(objectOnTimes)
        iOnTime = find(onTimes>objectOnTimes(j) & onTimes<objectOffTimes(j),1); % first saccade to this object
        if ~isempty(iOnTime)
            dwell{i}(j) = min(objectOffTimes(j),offTimes(iOnTime))-onTimes(iOnTime); % when the obj goes offscreen or the subj saccades away, whichever comes first
        else
            dwell{i}(j) = NaN;
        end
        isToTarget{i}(j) = objectIsTarget(j);
    end
    
    % Plot lines of objects seen in course of experiment
    if option_plot        
        for j=1:length(objectOnTimes)
            plot([objectOnTimes(j), objectOffTimes(j)]-y(i).eyelink.record_time,[i i],'r-');
        end      

        for j=1:length(onTimes)
            plot([onTimes(j) offTimes(j)]-y(i).eyelink.record_time,[i i]+0.1,'-');
        end
    end
    
    
end
% annotate plot
if option_plot
    ylim([0 numel(y)+1])
end
disp('Done!')