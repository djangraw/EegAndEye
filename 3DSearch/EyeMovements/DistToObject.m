function dist = DistToObject(x,positions,times)

% Calculate the distance between a given set of positions and the nearest
% visible objects at the given set of times.
%
% dist = DistToObject(x,positions,times)
%
% INPUTS:
% -x is a 3DSearch data struct
% -positions is an nx2 matrix in which the first column is x position and the
% second column is y position.
% -times is an n-element vector of the eyelink times at which locs should
% be evaluated.
% 
% OUTPUTS:
% -dist is an n-element column vector containing the distance from each
% point to the object visible at that time.  If no object was visible,
% dist=inf in that row.
%
% Created 8/16/11 by DJ.

% handle inputs
if length(positions)~=length(times)
    error('length of inputs ''positions'' and ''times'' must match!');
end

% Set up visible times info
visibletimes = GetVisibleTimes(x,'eyelink');
obj = visibletimes(:,1);
appeartimes = visibletimes(:,2);
disappeartimes = visibletimes(:,3);        
object_limits = x.eyelink.object_limits;

% set up dist
dist = Inf(length(times),1);

% Get distances
for i=1:numel(times)
    % Find all saccades made during each appearance and find distance from obj
    for j=1:numel(appeartimes)
        % Find the distance from each applicable saccade to the object's bounding box
        epoch_subset = find(times>appeartimes(j) & times < disappeartimes(j));
        t = times(epoch_subset); % times when these saccades were made
        for k=1:numel(epoch_subset)
            % Get the time when this object was last seen
            iLastLimit = find(object_limits(:,2)==obj(j) & object_limits(:,1)<t(k),1,'last');
            % Get object limits
            left = object_limits(iLastLimit,3);
            top = object_limits(iLastLimit,4);
            width = object_limits(iLastLimit,5);
            height = object_limits(iLastLimit,6);
            % Get saccade position
            px = positions(epoch_subset(k),1);
            py = positions(epoch_subset(k),2);
            % Get shortest distance from saccade position to object box
            if px<left
                dx = left-px;
            elseif px>left+width
                dx = px-(left+width);
            else 
                dx = 0;
            end
            if py<top
                dy = top-py;
            elseif py>top+height
                dy = py-(top+height);
            else
                dy = 0;
            end
            dist(epoch_subset(k)) = sqrt(dx*dx + dy*dy);
        end
    end
end