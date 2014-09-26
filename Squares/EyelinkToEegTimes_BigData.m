function eeg_times = EyelinkToEegTimes_BigData(eye_times,iSessions,BigData)

% eeg_times = EyelinkToEegTimes_BigData(eye_times,iSessions,BigData)
%
% - Converts Eyelink clock times to eeg clock times (both in # samples) 
% using linear interpolation. This usually results in a lag of no more than 
% 2ms (try running CheckSyncEvents(x) for details).
%
% INPUTS: 
% - BigData is a Squares data struct imported by import_squares_data and 
% combined by GetBigDataStruct.  All the necessary data can be extracted 
% from this struct.
% - eye_times is a vector of length n of the eyelink times you want to 
% convert to eeg times.
% - iSessions is a vector of length n of the session index of each
% eye_time.
% 
% OUTPUTS:
% - eeg_times is a vector of the same length as eye_times.
%
% Created 11/16/11 by DJ based on EyelinkToEegTimes.

all_iSessions = unique(iSessions);

eeg_times = nan(size(eye_times));
for i=1:numel(all_iSessions)
    % Get sync times for this session
    eye_syncs = BigData.sync.eyelink(BigData.sync.iSession==all_iSessions(i));
    eeg_syncs = BigData.sync.eeg(BigData.sync.iSession==all_iSessions(i));
    eeg_syncs = eeg_syncs + BigData.eeg.session_start_times(all_iSessions(i));
    % Error check
    if numel(eye_syncs) ~= numel(eeg_syncs)
        error('The number of sync signals sent by eyelink and received by eeg do not match!')
    end
    % Interpolate    
    isInSession = iSessions==all_iSessions(i);
    eeg_times(isInSession) = interp1(eye_syncs,eeg_syncs,eye_times(isInSession),'linear','extrap');    
end
