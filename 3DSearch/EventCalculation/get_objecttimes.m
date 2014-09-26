function object_times = get_objecttimes(object_limits)

% Generates a timestamp and event number for each time an object enters or
% exits the subject's view.
%
% object_times = get_objecttimes(object_limits)
% 
% INPUT object_limits is an nx6 matrix, where n is the number of frames in
% which an object was visible.  It is the output of get_objectlimits.m.
% OUTPUT object_times is an mx2 matrix, where m is the number of frames in
% which an object entered or exited the scene.
%
% Created 8/18/10 by DJ.
% Updated 8/19/10 by DJ - added Numbers.LAG_BETWEEN_FRAMES

% Setup
object_nums = unique(object_limits(:,2)); % all objects that were seen
times = object_limits(:,1); % extract times for easy access.
object_times = [];
GetNumbers;
threshold = Numbers.MAX_LAG_BETWEEN_FRAMES; % time (in ms) between sightings for something to be considered an exit and re-entrance.

% Get object events
for i = 1:numel(object_nums)
    ontimes = times(object_limits(:,2)==object_nums(i));
    object_times = [object_times; ontimes(1) Numbers.ENTERS + object_nums(i)];
    iGaps = find(diff(ontimes)>threshold);
    for j = 1:numel(iGaps)
        object_times = [object_times; ontimes(iGaps(j)) Numbers.EXITS + object_nums(i)];
        object_times = [object_times; ontimes(iGaps(j)+1) Numbers.ENTERS + object_nums(i)];
    end
    object_times = [object_times; ontimes(end) Numbers.EXITS + object_nums(i)];
end

% Sort into chronological order
[foo order] = sort(object_times(:,1));
object_times = object_times(order,:);