function ChangeSacToObjThresholds(subject,sessions,pixelThresholds,timeLimits,eventsRule,maxSeeTime,singleSuffix,comboSuffix)

% Changes the thresholds for detecting saccades to objects and recalculates
% EEGLAB events.
%
% ChangeSacToObjThresholds(subject,sessions,pixelThresholds,timeLimits,
% eventsRule,maxSeeTime,singleSuffix,comboSuffix)
%
% INPUTS:
% -subject is the subject number
% -sessions is a vector of all the sessions for that subject.
% -pixelThresholds is a vector of distances (in pixels, from the saccade 
% endpoint to the object) that a saccade must be within to be considered a
% saccade "to the object". The first saccade inside thresholds(1) will be 
% used. If that doesn't exist, the first saccade inside threshold 2 will be
% used (and so on). [Default: 50 100 150]
% -timeLimits is a 2-element vector with the start and end times, in ms,
% that are acceptable for a saccade to be counted.  [Default: 0 Inf]
%
% Created 2/23/11 by DJ.
% Updated 2/25/11 by DJ - added lots of inputs.  TO DO: update comments

% Handle defaults
if nargin<8 || isempty(comboSuffix)
    comboSuffix = '-all-filtered';
end
if nargin<7 || isempty(singleSuffix)
    singleSuffix = '-filtered';
end
if nargin<6 || isempty(maxSeeTime)
    maxSeeTime = 200; % max time after stim before an object is considered seen.
end
if nargin<5 || isempty(eventsRule)
    eventsRule = 'EarlySaccades';
end

% Add events
for i=1:numel(sessions)
    session = sessions(i);
    load(sprintf('3DS-%d-%d.mat',subject,session));
    x.eyelink.saccade_events = classify_saccades(x.eyelink.saccade_times,x.eyelink.saccade_positions,x.eyelink.object_limits,pixelThresholds,timeLimits);
    x = EyelinkToEegTimes(x);
    save(sprintf('3DS-%d-%d.mat',subject,session), 'x');
    
    AddEeglabEvents(subject,session,singleSuffix,eventsRule,maxSeeTime);
end
% Recombine sessions
CombineEeglabSessions(subject,sessions,singleSuffix,comboSuffix);