function [tPup, dPup] = RemoveSaccadePupilChanges(data,subject,sessions)

% Removes the pupil data changes that happen during each saccade.
%
% [tPup, dPup] = RemoveSaccadePupilChanges(data,subject,sessions)
%
% INPUTS: 
% - data is the h.data output from PlotEyeErps_MultiSession.m.
% - subject is the number of the subject.
% - sessions are the session numbers that have been imported using
% import_3DS_data_v3.
%
% OUTPUTS:
% - tPup is the target trial pupil response data with saccade changes 
% removed.  It should replace data.targetPupEpochs.
% - dPup is the distractor trial pupil response data with saccade changes 
% removed.  It should replace data.distractorPupEpochs.
%
% Created 11/4/11 by DJ.

for i=1:numel(sessions)
    load(sprintf('3DS-%d-%d',subject,sessions(i)))
    y(i) = x;
end

% Extract info
tPup = data.targetPupEpochs;
dPup = data.distractorPupEpochs;
times = data.epochTimes;
targetSessions = data.targetEventSessions;
targetOffsets = data.targetEventTimes;
distractorSessions = data.distractorEventSessions;
distractorOffsets = data.distractorEventTimes;

for i=1:size(tPup,1)
    tSacStart = y(sessions==targetSessions(i)).eyelink.saccade_start_times;
    tSacEnd = y(sessions==targetSessions(i)).eyelink.saccade_times;
    tEpoch = targetOffsets(i)+times;
    tMin = tEpoch(1);
    tMax = tEpoch(end);
    
    % Get saccade and event info
    tSacStart_inepoch = tSacStart(tSacStart > tMin & tSacStart < tMax);
    tSacEnd_inepoch = tSacEnd(tSacEnd > tMin & tSacEnd < tMax);
    if tSacEnd_inepoch(1) < tSacStart_inepoch(1)
        tSacStart_inepoch = [tMin; tSacStart_inepoch]; 
    end
    if tSacEnd_inepoch(end) < tSacStart_inepoch(end)
        tSacEnd_inepoch = [tSacEnd_inepoch; tMax];
    end
    
%     cla;
%     hold on;
%     plot(tEpoch,tPup(i,:),'b');    
%     ylabel('pupil size (a.u.)');
%     xlabel('time from object appearance (ms)')
    
    % Remove saccade changes
    for j=1:numel(tSacEnd_inepoch)
        iFirstIn = find(tEpoch>=tSacStart_inepoch(j),1);
        iFirstOut = find(tEpoch>=tSacEnd_inepoch(j),1);
        tPup(i,iFirstIn:iFirstOut-1) = tPup(i,iFirstIn); % flat during saccade
        tPup(i,iFirstOut:end) = tPup(i,iFirstOut:end)+tPup(i,iFirstIn)-tPup(i,iFirstOut); % remove discontinuity
    end
    
%     plot(tEpoch,tPup(i,:),'r');
%     PlotVerticalLines([tSacStart_inepoch],'g:');
%     PlotVerticalLines([tSacEnd_inepoch],'g-');
%     pause;
end

for i=1:size(dPup,1)
    tSacStart = y(distractorSessions(i)).eyelink.saccade_start_times;
    tSacEnd = y(distractorSessions(i)).eyelink.saccade_times;
    tEpoch = distractorOffsets(i)+times;
    tMin = tEpoch(1);
    tMax = tEpoch(end);
    
    % Get saccade and event info
    tSacStart_inepoch = tSacStart(tSacStart > tMin & tSacStart < tMax);
    tSacEnd_inepoch = tSacEnd(tSacEnd > tMin & tSacEnd < tMax);
    if tSacEnd_inepoch(1) < tSacStart_inepoch(1)
        tSacStart_inepoch = [tMin; tSacStart_inepoch]; 
    end
    if tSacEnd_inepoch(end) < tSacStart_inepoch(end)
        tSacEnd_inepoch = [tSacEnd_inepoch; tMax];
    end
    
%     cla;
%     hold on;
%     plot(tEpoch,dPup(i,:),'b');    
%     ylabel('pupil size (a.u.)');
%     xlabel('time from object appearance (ms)')
%     
    % Remove saccade changes
    for j=1:numel(tSacEnd_inepoch)
        iFirstIn = find(tEpoch>=tSacStart_inepoch(j),1);
        iFirstOut = find(tEpoch>=tSacEnd_inepoch(j),1);
        dPup(i,iFirstIn:iFirstOut-1) = dPup(i,iFirstIn); % flat during saccade
        dPup(i,iFirstOut:end) = dPup(i,iFirstOut:end)+dPup(i,iFirstIn)-dPup(i,iFirstOut); % remove discontinuity
    end
    
%     plot(tEpoch,dPup(i,:),'r');
%     PlotVerticalLines([tSacStart_inepoch],'g:');
%     PlotVerticalLines([tSacEnd_inepoch],'g-');
%     pause;
    
end