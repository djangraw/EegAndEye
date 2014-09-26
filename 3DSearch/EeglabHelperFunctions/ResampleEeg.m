function x = ResampleEeg(x,newsamplerate)

% -Takes the event/time fields in x.eeg and sets them to a new sample rate.
%
% x = ResampleEeg(x,newsamplerate)
% 
% -INPUT x is a 3DSearch data structure as imported by ImportData.
% -INPUT newsamplerate is the sampling rate in Hz to which you'd like to
% convert your EEG time data. We assume that the EEG data starts at time 0.
%
% Created 7/29/10 by DJ.
% Updated 8/19/10 by DJ - added saccade_events.
% Updated 8/31/10 by DJ - added record_time.
% Updated 12/16/10 by DJ - added blink_times.
% Updated 12/14/11 by DJ - added sync.eeg.

% find the factor by which you should multiply eeg time data
oldsamplerate = x.eeg.eventsamplerate;
ratio = newsamplerate/oldsamplerate;

% Get fields
nFields = numel(fieldnames(x.eeg))-1; % exclude eventsamplerate as a field
if isfield(x,'sync')
    nFields = nFields+1; 
end
nResampled = 0;

% Resample each field, if it exists
if isfield(x.eeg,'sync_events')
    x.eeg.sync_events(:,1) = x.eeg.sync_events(:,1) * ratio;
    nResampled = nResampled+1;
end
if isfield(x.eeg,'object_events')
    x.eeg.object_events(:,1) = x.eeg.object_events(:,1) * ratio;
    nResampled = nResampled+1;
end
if isfield(x.eeg,'saccade_events')
    x.eeg.saccade_events(:,1) = x.eeg.saccade_events(:,1) * ratio;
    nResampled = nResampled+1;
end
if isfield(x.eeg,'object_lifetimes')
    x.eeg.object_lifetimes = x.eeg.object_lifetimes * ratio;
    nResampled = nResampled+1;
end
if isfield(x.eeg,'button_times')
    x.eeg.button_times = x.eeg.button_times * ratio;
    nResampled = nResampled+1;
end
if isfield(x.eeg,'blink_times')
    x.eeg.blink_times = x.eeg.blink_times * ratio;
    nResampled = nResampled+1;
end
if isfield(x.eeg,'record_time')
    x.eeg.record_time = x.eeg.record_time * ratio;
    nResampled = nResampled+1;
end
if isfield(x.eeg,'load_time')
    x.eeg.load_time = x.eeg.load_time * ratio;
    nResampled = nResampled+1;
end
if isfield(x.eeg,'start_time')
    x.eeg.start_time = x.eeg.start_time * ratio;
    nResampled = nResampled+1;
end
if isfield(x.eeg,'end_time')
    x.eeg.end_time = x.eeg.end_time * ratio;
    nResampled = nResampled+1;
end
if isfield(x,'sync')
    x.sync.eeg = x.sync.eeg * ratio;
    nResampled = nResampled+1;
end

% set new samplerate field
x.eeg.eventsamplerate = newsamplerate;

% display results
fprintf('Resampled %d of %d fields in x.eeg\n',nResampled,nFields);

