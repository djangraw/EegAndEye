function delay = CheckSyncEvents(x)

% delay = CheckSyncEvents(x);
%
% - Extracts and compares the sync event times according to the eye tracker
% and the EEG to check for consistency.
% - Sync events are sent every couple of seconds from the eyelink (which
% records the time each was sent) to the EEG (which records the time it was
% received). 
% - Sync events will be used to translate eyelink times into EEG times for
% analysis (e.g. by function EyelinkToEegTimes).
%
% Created 6/14/10 by DJ (as part of CheckData).
% Updated 7/28/10 by DJ (made into its own program).
% Updated 7/29/10 by DJ (changed events field back to eyelink).
% Updated 12/6/13 by DJ (added delay output, cleaned up).

% Extract info from the struct we made (for easier access)
eye = x.eyelink.sync_events; % should start from the 'START_RECORDING' broadcast
eeg = x.eeg.sync_events; % should start from the 'START_RECORDING' broadcast

% Take care of sampling rates, offsets
eeg(:,1) = eeg(:,1)*1000/x.eeg.eventsamplerate; % to get time in ms instead of samples (EEG sampling rate = 2048)
eeg(:,1) = eeg(:,1) - eeg(1,1); % make time of first event t=0
eye(:,1) = eye(:,1) - eye(1,1); % make time of first event t=0

% Make sure events are the same
if size(eeg,1) ~= size(eye,1) || sum(eeg(:,2)~=eye(:,2)) > 0 % number of events that aren't in exactly the right spot
    error('Sync events don''t match up!');
else
    disp('All event types match. Checking timing...')
    delay = diff(eeg(:,1))-diff(eye(:,1)); % subtract the time between timestamps
    subplot(2,1,1)
    hist(delay);
    xlabel('Delay (Time Received - Time Sent) in ms')
    ylabel('# events')
    title('Event Timing Consistency Check');
    subplot(2,1,2);
    plot(eeg(2:end,1)/1000,delay, '.');
    xlabel('Event time (s)');
    ylabel('Delay in ms');
    fprintf('average absolute value of delay is %g\n', mean(abs(delay)));
end
