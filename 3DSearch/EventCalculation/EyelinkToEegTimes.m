function output = EyelinkToEegTimes(varargin)

% x = EyelinkToEegTimes(x)
% eeg_times = EyelinkToEegTimes(eye_times,x)
% eeg_times = EyelinkToEegTimes(eye_times,eye_syncs,eeg_syncs)
%
% - Converts Eyelink clock times to eeg clock times (both in # samples) 
% using linear interpolation. This usually results in a lag of no more than 
% 2ms (try running CheckSyncEvents(x) for details).
%
% INPUTS: 
% - x is a 3DSearch data struct imported by ImportData.  All the
% necessary data can be extracted from this struct.
% - eye_times is a vector of the eyelink times you want to convert to 
% eeg times.
% - eye_syncs and eeg_syncs are the times when sync signals were
% sent to and received by the EEG (respectively).  They must have the same 
% number of elements.
% 
% OUTPUTS:
% - eeg_times is a vector of the same length as eye_times.
%
% Format 1 will convert all the events in x.eyelink to EEG times and place
% them in x.eeg.
% Format 2 is just like format 3, but will extract all the other info
% automatically from the struct x.
%
% Created 7/27/10 by DJ.
% Updated 7/28/10 by DJ (debugged)
% Updated 7/29/10 by DJ (added input options)
% Updated 8/3/10 by DJ (switched to linear interpolation)
% Updated 8/19/10 by DJ (added saccade_events)
% Updated 8/31/10 by DJ (added record_time)
% Updated 10/18/10 by DJ (added blink_times, comments)
% Updated 7/28/11 by DJ (added brake_button_times)
% Updated 8/10/11 by DJ (added check that brake_button_times field exists)
% Updated 11/11/11 by DJ (removed eyetimes)
% Updated 11/16/11 by DJ (added squares struct compatability (sync field))
% Updated 12/5/13 by DJ (added leader_slow_times)

if nargin==1 % SPECIAL CASE - input format x = EyelinkToEegTimes(x)
    x = varargin{1};
    
    % Call this function recursively to update each field
    x.eeg.record_time = EyelinkToEegTimes(x.eyelink.record_time,x);
    x.eeg.load_time = EyelinkToEegTimes(x.eyelink.load_time,x);
    x.eeg.start_time = EyelinkToEegTimes(x.eyelink.start_time,x);
    x.eeg.end_time = EyelinkToEegTimes(x.eyelink.end_time,x);

    x.eeg.object_lifetimes(:,1) = EyelinkToEegTimes(x.eyelink.object_lifetimes(:,1),x);
    x.eeg.object_lifetimes(:,2) = EyelinkToEegTimes(x.eyelink.object_lifetimes(:,2),x);
    x.eeg.object_events = x.eyelink.object_events; % column 2 (event codes) will be the same
    x.eeg.object_events(:,1) = EyelinkToEegTimes(x.eyelink.object_events(:,1),x);
    x.eeg.saccade_events = x.eyelink.saccade_events; % column 2 (event codes) will be the same
    x.eeg.saccade_events(:,1) = EyelinkToEegTimes(x.eyelink.saccade_events(:,1),x);
    
    x.eeg.button_times = EyelinkToEegTimes(x.eyelink.button_times,x);
    if isfield(x.eyelink,'brake_button_times')
        x.eeg.brake_button_times = EyelinkToEegTimes(x.eyelink.brake_button_times,x);
    else
        x.eyelink.brake_button_times = [];
        x.eeg.brake_button_times = [];        
    end
    if isfield(x.eyelink,'leader_slow_times')
        x.eeg.leader_slow_times = EyelinkToEegTimes(x.eyelink.leader_slow_times,x);
    else
        x.eyelink.leader_slow_times = [];
        x.eeg.leader_slow_times = [];        
    end
    x.eeg.blink_times = EyelinkToEegTimes(x.eyelink.blink_times,x);
    
    output = x;
    
else % NORMAL CASE - just convert the given times
    % Handle Inputs    
    if nargin==2 % if the x input option was used: eeg_times = EyelinkToEegTimes(eye_times,x)
        % extract the requested input data from x
        eye_times = varargin{1};
        x = varargin{2};
        if isfield(x,'eyelink'); % 3DSearch file
            eye_syncs = x.eyelink.sync_events(:,1);
            eye_num = x.eyelink.sync_events(:,2);
            eeg_syncs = x.eeg.sync_events(:,1);
            eeg_num = x.eeg.sync_events(:,2);
        elseif isfield(x,'sync'); % squares file
            eye_syncs = x.sync.eyelink;
            eeg_syncs = x.sync.eeg;
        end
    elseif nargin>=3
        eye_times = varargin{1};
        eye_syncs = varargin{2};
        eeg_syncs = varargin{3};
    end

    % eye_syncs and eeg_syncs should be vectors, but if they are the
    % sync_events (nx2) matrices, convert them!
    if size(eye_syncs,2)==2
        eye_syncs=eye_syncs(:,1);
    end
    if size(eeg_syncs,2)==2
        eeg_syncs=eeg_syncs(:,1);
    end
           
    
    % Error check
    if numel(eye_syncs) ~= numel(eeg_syncs)
        error('The number of sync signals sent by eyelink and received by eeg do not match!')
    end

    % Linearly interpolate new times
    eeg_times = interp1(eye_syncs,eeg_syncs,double(eye_times),'linear','extrap');
    % Send to output
    output = eeg_times;
end