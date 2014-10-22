function EEG = NEDE_RejectSession(EEG,iSessions)

% EEG = NEDE_RejectSession(EEG,iSessions)
%
% INPUTS:
% - EEG is a continuous (not epoched) eeglab data struct
% - iSessions is a vector of indices of the sessions you want to reject.
%
% OUTPUTS:
% - EEG is the input struct, but with the offending sessions removed.
%
% Created 10/22/14 by DJ.

% get boundary events
types = {EEG.event.type};
isBound = strcmp(types,'boundary');
% get boundary event latencies (in samples)
boundLat = [1, EEG.event(isBound).latency, EEG.pnts];

% get start and end time for each session
regions = nan(numel(iSessions),2);
for i=1:numel(iSessions)    
    regions(i,1) = boundLat(iSessions(i));
    regions(i,2) = boundLat(iSessions(i)+1);
end

% reject times
EEG = eeg_eegrej( EEG, regions);
