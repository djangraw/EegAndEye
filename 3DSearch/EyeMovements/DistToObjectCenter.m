function dist = DistToObjectCenter(x,positions,times,uselastvisiblespot)

% Calculate the distance between a given set of positions and the nearest
% visible objects at the given set of times.
%
% dist = DistToObjectCenter(x,positions,times,uselastvisiblespot)
%
% INPUTS:
% -x is a 3DSearch data struct
% -positions is an nx2 matrix in which the first column is x position and the
% second column is y position.
% -times is an n-element vector of the eyelink times at which locs should
% be evaluated.
% -uselastvisiblespot is a binary value indicating whether, if an object
% isn't on screen, the last place an obect was visible should be used.
% 
% OUTPUTS:
% -dist is an n-element column vector containing the distance from each
% point to the center of the object visible at that time.  If no object was 
% visible: if uselastvisiblespot, dist=inf, otherwise it's the distance to
% the last spot where an object was visible.
%
% Created 8/26/13 by DJ based on DistToObj.m.

if nargin<4 || isempty(uselastvisiblespot)
    uselastvisiblespot = 1;
end

% handle inputs
if length(positions)~=length(times)
    error('length of inputs ''positions'' and ''times'' must match!');
end

% Set up visible times info
visibletimes = GetVisibleTimes(x,'eyelink');
obj = visibletimes(:,1);
appeartimes = visibletimes(:,2);
disappeartimes = visibletimes(:,3);        
% kluge for uselastvisiblespot (pretend the object never disappeared!)
if uselastvisiblespot
    disappeartimes(1:end-1) = appeartimes(2:end);
    disappeartimes(end) = inf;
end

object_limits = x.eyelink.object_limits;

% set up dist
dist = Inf(length(times),1);

% ---Get distances
% Find all saccades made during each appearance and find distance from obj
for j=1:numel(appeartimes)
    % Find the distance from each applicable saccade to the object's bounding box
    epoch_subset = find(times>appeartimes(j) & times <= disappeartimes(j));
    t = times(epoch_subset); % times when these saccades were made
    for k=1:numel(epoch_subset)
        % Get the time when this object was last seen
        iLastLimit = find(object_limits(:,2)==obj(j) & object_limits(:,1)<t(k),1,'last');
        % Get object limits
        left = object_limits(iLastLimit,3);
        top = object_limits(iLastLimit,4);
        width = object_limits(iLastLimit,5);
        height = object_limits(iLastLimit,6);
        % Get object center
        cx = left + width/2;
        cy = top + height/2;
        % Get saccade position
        px = positions(epoch_subset(k),1);
        py = positions(epoch_subset(k),2);
        % Get distance from saccade position to object center
        dx = cx - px;
        dy = cy - py;
        dist(epoch_subset(k)) = sqrt(dx*dx + dy*dy);
    end
end

% Fix any early infinity values (KLUSGE! IS THIS THE RIGHT THING TO DO?)
early_subset = find(isinf(dist));
if ~isempty(early_subset) && uselastvisiblespot
    
    % Get the time when the first object was first seen
    [~,iFirstLimit] = min(object_limits(:,1));
    % Get object limits
    left = object_limits(iFirstLimit,3);
    top = object_limits(iFirstLimit,4);
    width = object_limits(iFirstLimit,5);
    height = object_limits(iFirstLimit,6);
    % Get object center
    cx = left + width/2;
    cy = top + height/2;
    
    for k=1:numel(early_subset)
        % Get saccade position
        px = positions(early_subset(k),1);
        py = positions(early_subset(k),2);
        % Get distance from saccade position to object center
        dx = cx - px;
        dy = cy - py;
        dist(early_subset(k)) = sqrt(dx*dx + dy*dy);
    end
end
