function TEMP_AddRecordingStartTime(subject,sessions)

% Update the .mat files with the time when eyelink recording started.
%
% Created 7/2/12 by DJ for one-time use.

for i=1:numel(sessions)
    filename = sprintf('sq-%d-%d',subject,sessions(i));
    fprintf('%d...',i);
    load(filename);
    x.recording_start_time = find_first_word(x.eyeFilename,'START','%d');
    save(filename,'x');
end
fprintf('all...');
saveBehaviorData(subject,sessions);
disp('Done!')
