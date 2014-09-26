function BigData = GetBigDataStruct(subject,sessions,experimentType)

% Loads all the data for a single subject and puts it in one big struct
%
% Created 11/16/11 by DJ.

% Load
switch experimentType
    case '3DSearch'
        for i=1:numel(sessions)
            load(sprintf('3DS-%d-%d.mat',subject,sessions(i)));
            y(i) = x;
        end
    case 'Squares'
        for i=1:numel(sessions)
            load(sprintf('sq-%d-%d.mat',subject,sessions(i)));
            y(i) = x;
        end
end

% Combine
BigData = combine_structs(y,[],experimentType);

% Add EEG information
eegfilename = sprintf('sq-%d-all-filtered.set',subject);
EEG = pop_loadset(eegfilename);

iBoundary = strmatch('boundary',{EEG.event(:).type});
session_start_times = [0 EEG.event(iBoundary).latency];
BigData.eeg.combined_filename = eegfilename;
BigData.eeg.session_start_times = session_start_times;
