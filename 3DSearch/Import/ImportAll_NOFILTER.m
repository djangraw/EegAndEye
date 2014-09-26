% Get info from GUI fields

handles.subject = 19;
handles.sessions = 2:11;
handles.eegsessions = 1:10;
handles.filetype = 'Sensorium-2011';
handles.eventsrule = 'OddballTask';
handles.pixelthresh = 200;
handles.timelimits = [0 Inf];
handles.maxseetime = 2000;
handles.singlesuffix = '-raw';
handles.combosuffix = '-all-raw';

for i=1:numel(handles.sessions)
    % Run functions as in ImportData
    Import_3DS_Data_v3(handles.subject,handles.sessions(i),handles.eegsessions(i),handles.filetype,handles.pixelthresh,handles.timelimits);
    ImportToEeglab(handles.subject,handles.sessions(i),handles.eegsessions(i),handles.filetype, '3DSearch',handles.singlesuffix,[NaN NaN],[NaN NaN],NaN,true); % no filters and include EOG
    AddEeglabEvents(handles.subject,handles.sessions(i),handles.singlesuffix,handles.eventsrule,handles.maxseetime);
end
%Combine sessions
CombineEeglabSessions(handles.subject,handles.sessions,handles.singlesuffix,handles.combosuffix);
