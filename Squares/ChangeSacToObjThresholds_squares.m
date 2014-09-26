function ChangeSacToObjThresholds_squares(subject,sessions,pixelThreshold,eventsRule,singleSuffix,comboSuffix)

% Changes the thresholds for detecting saccades to objects and recalculates
% EEGLAB events.
%
% ChangeSacToObjThresholds_squares(subject,sessions,pixelThreshold,eventsRu
% le,singleSuffix,comboSuffix)
%
% INPUTS:
% -subject is the subject number
% -sessions is a vector of all the sessions for that subject.
% -pixelThreshold is a scalar indicating the distance(in pixels, from the  
% saccade endpoint to the square center) that a saccade must be within to 
% be considered a saccade "to the square".
%
% Created 10/31/11 by DJ based on program ChangeSacToObjThresholds_squares.

% Handle defaults
if nargin<6 || isempty(comboSuffix)
    comboSuffix = '-all-filtered';
end
if nargin<5 || isempty(singleSuffix)
    singleSuffix = '-filtered';
end
if nargin<4 || isempty(eventsRule)
    eventsRule = 'Squares';
end

Constants = GetSquaresConstants;

% Add events
for i=1:numel(sessions)
    session = sessions(i);
    filename = sprintf('sq-%d-%d.mat',subject,session);
    fprintf('Updating saccades in file %s...\n',filename);
    load(filename);
    % change pixel threshold
    x.pixel_threshold = pixelThreshold;
    x.eyelink.saccade_events = classify_squares_saccades(x);        
    
    % Update reaction times - from import_squares_data.m
    median_iti = median(diff(x.trial.start_time));
    completion_time = zeros(size(x.trial.start_time));
    response_time = zeros(size(x.trial.start_time));
    reaction_time = zeros(size(x.trial.start_time));
    is_target_response = nan(size(x.trial.start_time));
    for i=1:numel(x.trial.start_time)
        if x.trial.is_target_trial(i)
            ct = x.saccade.end_time(find(x.saccade.end_time>x.trial.start_time(i) & x.saccade.class==Constants.COMPLETION, 1));
        elseif ~x.trial.is_right_cross(i)
            ct = x.saccade.end_time(find(x.saccade.end_time>x.trial.start_time(i) & x.saccade.squarenum==numel(Constants.SQUARE_X), 1));
        else
            ct = x.saccade.end_time(find(x.saccade.end_time>x.trial.start_time(i) & x.saccade.squarenum==1, 1));
        end
        if ~isempty(ct)
            completion_time(i) = ct;
        end

        iFirstButton = find(x.button.time>x.trial.start_time(i),1);
        if isempty(iFirstButton) || x.button.time(iFirstButton)-x.trial.start_time(i)>median_iti
            response_time(i) = NaN;
            reaction_time(i) = NaN;
        else
            response_time(i) = x.button.time(iFirstButton);
            reaction_time(i) = response_time(i)-completion_time(i);
            is_target_response(i) = (x.button.number(iFirstButton)==Constants.TARGET_BUTTON);
        end    
    end
    x.trial.completion_time = completion_time;
    x.trial.response_time = response_time;
    x.trial.reaction_time = reaction_time;
    
    % save results
    save(sprintf('sq-%d-%d.mat',subject,session), 'x');
    % add to eeglab file
    AddEeglabEvents(subject,session,singleSuffix,eventsRule);
end
% Recombine sessions
CombineEeglabSessions(subject,sessions,singleSuffix,comboSuffix,'Squares');