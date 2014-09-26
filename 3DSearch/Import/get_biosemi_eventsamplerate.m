function sampling_rate = get_biosemi_eventsamplerate(filename)

% sampling_rate = get_biosemi_eventsamplerate(filename)
%
% Given a biosemi filename (.bdf) in the current path, loads the info about
% the file and uses it to find the sampling rate of the events channel (in
% Hz).
%
% Created 6/16/10 by DJ.

% Open file
dat = sopen(filename);

% extract sampling rate for event channel (we assume it's the last channel)
sampling_rate = dat.AS.SampleRate(end); 

% Clean up
sclose(dat);