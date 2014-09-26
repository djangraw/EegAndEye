function ts_events = GetRecognizableTimes(object_limits,sac_time, saccade_positions,cutoff_frac, cutoff_size,cutoff_ecc)

% Finds the times at which objects first became recognizable, as defined by certain cutoffs.
%
% ts_events = GetRecognizableTimes(object_limits,sac_time, ...
% saccade_positions,cutoff_frac, cutoff_size,cutoff_ecc)
% 
% INPUTS:
% - object_limits is an nx7 matrix, where n is the number of frames in
% which an object was visible.  Each row is [timestamp, object, left, top,
% width, height, fractionvisible].
% - saccade_times and saccade_positions are the timestamps and x,y
% positions of every saccade made during the session.
% - cutoff_frac, cutoff_size, and cutoff_ecc are the minimum fraction 
% visible (as reported by a Unity replay), minimum size (defined as the 
% greater of width,length), and maximum eccentricity (defined as the
% shortest dist in pixels between the eyes and the object), respectively.
% 
% OUTPUTS:
% - ts_events is an mx2 matrix, where m is the number of "became 
% recognizable" events found. Each row is [timestamp eventnumber].  Event
% numbers are explained in GetNumbers.m.
%
% Created 8/25/10 by DJ.

% Parse out inputs
lim_time = object_limits(:,1); % time at which object limit was recorded
lim_obj =  object_limits(:,2); % object number
lim_lims =  object_limits(:,3:6); % left, top, width, height
lim_frac =  object_limits(:,7); % fraction of object visible

sac_pos = saccade_positions(:,2:3);

nLim = numel(lim_time);
nObj = max(lim_obj);
recognizable = zeros(nLim,nObj);

% Main loop
for i=1:nLim
    % Find object info
    t = lim_time(i);    
    obj_cx = lim_lims(i,1) + lim_lims(i,3)/2;
    obj_cy = lim_lims(i,2) + lim_lims(i,4)/2;
    % Find saccade info
    last_sac = find(sac_time <= t,1,'last');
    x = sac_pos(last_sac,1);
    y = sac_pos(last_sac,2);
    % Find size, eccentricity and frac visible
    obj_frac = lim_frac(i);
    obj_size = max(lim_lims(i,3:4));
    obj_ecc = sqrt((x-obj_cx)^2 + (y-obj_cy)^2);
    % Compare to cutoffs
    if obj_frac >= cutoff_frac && obj_size >= cutoff_size && obj_ecc <= cutoff_ecc
        recognizable(i,lim_obj(i)) = 1;
    end
end

% Find times when an object was recognizable for >x ms after not being
% recognizable for y ms before that.
ts_events = [];
x = 100;
y = 100;
GetNumbers;
for j = 1:nObj
    change = diff(recognizable(:,j));
    uptimes = lim_time(find(change==1));
    downtimes = lim_times(find(change==-1));
    if numel(uptimes)>numel(downtimes)
        downtimes = [downtimes; lim_time(end)];
    elseif numel(downtimes)>numel(uptimes)
        uptimes = [uptimes; lim_time(1)];
    end
    
    if downtimes(1)-uptimes(1) > x
        ts_events = [ts_events; uptimes(1) Numbers.ENTERS+j];
    end
    for k=2:numel(uptimes)
        if downtimes(k)-uptimes(k)>x && uptimes(k)-downtimes(k-1)>y
            [ts_events; uptimes(k) Numbers.ENTERS+j];
        end
    end
    
end
