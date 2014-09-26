function [sacstart_urevent, sacend_urevent, urevent_sacnum] = MatchUpSquaresEvents(subject,sessions)

% MatchUpEvents(subject,sessions)
%
% [sacstart_urevent, sacend_urevent, urevent_sacnum] = MatchUpSquaresEvents(subject,sessions)
%
% INPUTS:
% - subject is a scalar indicating the subject number in the filenames.
% - sessions is an n-element vector indicating the session numbers of the 
% datafiles to be analyzed.
%
% OUTPUTS:
% - sacstart_urevent is an n-element vector of cells, in which each cell is
% an m-element vector of urevent indices.  sacstart_urevent(i) is the index
% of the urevent corresponding to x.saccade.start_time(i).
% - sacend_urevent is an n-element vector of cells, in which each cell is
% an m-element vector of urevent indices.  sacend_urevent(i) is the index
% of the urevent corresponding to x.saccade.end_time(i).
% - urevent_sacnum is an n-element vector of cells, in which each cell is a 
% p-element vector of x.saccade indices.  urevent_sacnum(j) is the index of
% the saccade corresponding to EEG.urevent(j).
% 
% Created 10/31/11 by DJ.
% Updated 11/16/11 by DJ - to match new data structures' sync field

% startup
prefix = 'sq';
eeg_suffix = '-filtered';
Constants = GetSquaresConstants;
sacstart_urevent = cell(size(sessions));
sacend_urevent = cell(size(sessions));
urevent_sacnum = cell(size(sessions));

% Main loop
for i=1:numel(sessions)
    % set up
    eegfilename = sprintf('%s-%d-%d%s.set',prefix,subject,sessions(i),eeg_suffix);
    matfilename = sprintf('%s-%d-%d.mat',prefix,subject,sessions(i));
    fprintf('Matching up events for files %s and %s...\n', matfilename, eegfilename)
    % Load data
    EEG = pop_loadset('filename',eegfilename);    
    load(matfilename); % load trial structure x
    % Convert times
    sacstart_eegtime = EyelinkToEegTimes(x.saccade.start_time,x.sync.eyelink,x.sync.eeg)/x.eeg.eventsamplerate;
    sacend_eegtime = EyelinkToEegTimes(x.saccade.end_time,x.sync.eyelink,x.sync.eeg)/x.eeg.eventsamplerate; 
    class = x.saccade.class;
    class(isnan(class)) = Constants.OTHER;
    % Find matching events
    sacstart_urevent{i} = nan(size(sacstart_eegtime));
    sacend_urevent{i} = nan(size(sacend_eegtime));
    urevent_sacnum{i} = nan(size(EEG.urevent));    
    for j=1:numel(sacstart_eegtime)
        ssu = find([EEG.urevent(:).init_time]==sacstart_eegtime(j) & ...
            [EEG.urevent(:).type]==Constants.SACCADESTART_BASE+class(j), 1);
        seu = find([EEG.urevent(:).init_time]==sacend_eegtime(j) & ...
            [EEG.urevent(:).type]==Constants.SACCADEEND_BASE+class(j), 1);        
        if ~isempty(ssu) && ~isempty(seu)
            sacstart_urevent{i}(j) = ssu;    
            sacend_urevent{i}(j) = seu;
            urevent_sacnum{i}([ssu seu]) = j;
        end
    end
    
end
