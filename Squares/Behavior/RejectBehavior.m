function [isBadTrial, isBadSaccade] = RejectBehavior(x,rule)

% [isBadTrial, isBadSaccade] = RejectBehavior(x,rule)
%
% INPUTS:
% - x is a squares data file imported using import_squares_data.
% - rule is a string or cell array of strings indicating which type(s) of 
% trials you want to mark for rejection. Options: {'skip','skipped_ends',
% 'backward','errant','blink','early_button','late_button','early_saccade'}
% (see code for details).
%
% OUTPUTS:
% - isBadTrial is a vector of binary values, one per trial.  If
% isBadTrial(i) is true, trial i should be rejected based on the specified
% metric.
% - isBadTrial is a vector of binary values, one per saccade.  If
% isBadTrial(i) is true, saccade i should be rejected based on the 
% specified metric.
%
% Created 4/30/12 by DJ.
% Updated 5/8/12 by DJ - added early_button and late_button options
% Updated 6/18/12 by DJ - allow blank input rule
% Updated 8/2/12 by DJ - allow cell input, added early_saccade rule

% Declare defaults
if nargin<2
    rule = 'skip';
end

% Handle cell input by calling recursively
if iscell(rule)
    isBT = zeros(numel(x.trial.start_time), numel(rule));
    isBS = zeros(numel(x.saccade.start_time), numel(rule));
    for i=1:numel(rule)
        [isBT(:,i), isBS(:,i)] = RejectBehavior(x,rule{i});
    end
    isBadTrial = any(isBT,2);
    isBadSaccade = any(isBS,2);
    return;
end

% Set up
nTrials = numel(x.trial.start_time);
nSaccades = numel(x.saccade.start_time);
nSquares = size(x.trial.is_target_color,2);
isBadTrial = false(nTrials,1);
isBadSaccade = false(nSaccades,1);
Constants = GetSquaresConstants;

% return now if no rule is requested
if isempty(rule)
    return;
end

% Find trials and saccades to mark for rejection
for i=1:nTrials
    theseSaccades = find(x.saccade.trialnum==i);
    switch rule
        case {'backward_saccade' 'backward'}
            if any(x.saccade.class(theseSaccades) == Constants.BACKWARD) % if the subject looked backward in the sequence at any point
                isBadTrial(i) = true;
                isBadSaccade(theseSaccades) = true;
            end
        case {'skipped_square' 'skip'}
            if ~isempty(setdiff(1:nSquares,x.saccade.squarenum(theseSaccades))) % if any of the squares were not viewed
                isBadTrial(i) = true;
                isBadSaccade(theseSaccades) = true;
            end            
        case 'skipped_ends'
            if ~isempty(setdiff([0 nSquares+1],x.saccade.squarenum(theseSaccades))) % if the cross or circle were not viewed
                isBadTrial(i) = true;
                isBadSaccade(theseSaccades) = true;
            end
        case {'errant_saccade' 'errant'} % if any saccade did not land on a square (very strict!!!)
            if any(isnan(x.saccade.squarenum(theseSaccades)))
                isBadTrial(i) = true;
                isBadSaccade(theseSaccades) = true;
            end    
        case 'blink' % if a blink took place at any point during a trial
            if any(x.blink.end_time>x.trial.fix_time(i) & x.blink.start_time <x.trial.circle_time(i))
                isBadTrial(i) = true;
                isBadSaccade(theseSaccades) = true;
            end 
        case 'early_button' 
            if x.trial.reaction_time(i)<0 % if the subject responded before the trial ended
                isBadTrial(i) = true;
                isBadSaccade(theseSaccades) = true;
            end 
        case 'late_button' 
            if i<nTrials && x.trial.response_time(i)>x.trial.fix_time(i+1) % if subject responded after the next trial started
                isBadTrial(i) = true;
                isBadSaccade(theseSaccades) = true;
            end 
        case 'wrong_button'
            if ~x.trial.is_correct_response(i)
                isBadTrial(i) = true;
                isBadSaccade(theseSaccades) = true;
            end
        case 'early_saccade'
            if x.trial.is_right_cross(i)
                xdist = abs(x.saccade.end_position(theseSaccades,1)-Constants.SQUARE_X(end));
                ydist = abs(x.saccade.end_position(theseSaccades,2)-Constants.SQUARE_Y(end));                
            else
                xdist = abs(x.saccade.end_position(theseSaccades,1)-Constants.SQUARE_X(1));
                ydist = abs(x.saccade.end_position(theseSaccades,2)-Constants.SQUARE_Y(1));
            end
            dist = sqrt(xdist.^2 + ydist.^2);
            if any(dist<x.pixel_threshold & x.saccade.end_time(theseSaccades)<x.trial.start_time(i))
                isBadTrial(i) = true;
                isBadSaccade(theseSaccades) = true;
            end
        otherwise
            error('Rule ''%s'' not recognized!',rule)
    end
end