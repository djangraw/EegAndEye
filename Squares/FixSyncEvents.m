function [sync,eeg_events] = FixSyncEvents(x,eeg_events_in)

% Hard-coded sync event fixes for improperly recorded sessions.  
% 
% [sync,eeg_events] = FixSyncEvents(x,eeg_events_in)
%
% INPUTS:
% - x is a squares data struct
% - eeg_events_in is an nx2 matrix where the first column contains the
% times at which sync events were recorded in the eeg data and the second
% column contains the corresponding event numbers.
% 
% OUTPUTS:
% - sync is a struct with fields 'eyelink' and 'events' containing the
% times at which sync events were recorded in the eyelink data and the
% corresponding event numbers.
% - eeg_events is a fixed version of eeg_events_in.
%
% Created 2/13/12 by DJ.
% Updated 4/12/13 by DJ - added sqfix check

% Hard-code in individual fixes
if x.subject==8 && x.session==1 && strcmp(x.eegFilename(1:2),'sq')
    fprintf('HARD-CODED SYNC FIX FOR SUBJ %d, SESS %d',x.subject,x.session)
    sync.events = x.sync.events(3:end);
    sync.eyelink = x.sync.eyelink(3:end);
    eeg_events = eeg_events_in;
else
    sync = x.sync;
    eeg_events = eeg_events_in;
end
