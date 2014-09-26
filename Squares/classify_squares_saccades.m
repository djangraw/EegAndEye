function [saccade_trialnum, saccade_squarenum, saccade_class, saccade_distance] = classify_squares_saccades(sq,times,positions)

% [saccade_trialnum, saccade_squarenum, saccade_class, saccade_distance] = classify_squares_saccades(sq)
% ... = classify_squares_saccades(sq,times,positions)
% INPUTS:
% -sq is a squares data struct imported using import_squares_data.
% -times is an n-element vector of (eyelink) times at which saccades began.
% [default: sq.saccade.start_time]
% -positions is a size nx2 matrix of (screen) positions at which saccades
%  ended. [default: sq.saccade.end_position]
%
% OUTPUTS:
% -saccade_trialnum is a record of the trial during which this saccade took
%  place. It is an n-element vector, where n is the number of saccades.
% -saccade_squarenum is a record of which square number (increasing 
%  right-to-left) each saccade ended on. It is an n-element vector, where n
%  is the number of saccades.
% -saccade_class is a record of which type of saccade was made based on the
%  identity of the square and the classes of past saccades in this trial.
%  It is an n-element vector, where n is the number of saccades.  See
%  GetSquaresConstants for an explanation of each saccade class.
% -saccade_distance is a record of 
%
% Created 10/17/11 by DJ.
% Updated 10/19/11 by DJ - fixed WITHIN and BACKWARDS bugs, added output 
%  saccade_trialnum, restored dy
% Updated 11/28/11 by DJ - added optional inputs so program will work with
%  fixations too
% Updated 1/4/12 by DJ - added <= and recognizes saccades during endcircle
%  time.

if nargin<2 || isempty(times)
    times = sq.saccade.start_time;
end
if nargin<3 || isempty(positions)
    positions = sq.saccade.end_position;
end
    

Constants = GetSquaresConstants();
[saccade_squarenum,saccade_class,saccade_trialnum,saccade_distance] = ...
    deal(nan(size(times)));

for i=1:numel(sq.trial.start_time)
    % Check for saccades to fixation cross
    fixation_saccades = find(times>sq.trial.fix_time(i) & times<=sq.trial.start_time(i));
    saccade_trialnum(fixation_saccades) = i;
    for j=1:numel(fixation_saccades)
        iSac = fixation_saccades(j);
        if sq.trial.is_right_cross(i)
            dx = Constants.RIGHTCROSS_X-positions(iSac,1);
            dy = Constants.RIGHTCROSS_Y-positions(iSac,2);
        else
            dx = Constants.LEFTCROSS_X-positions(iSac,1);
            dy = Constants.LEFTCROSS_Y-positions(iSac,2);
        end
        disto = sqrt(dx*dx + dy*dy);
        if disto<sq.pixel_threshold
            saccade_distance(iSac) = disto;
            if sq.trial.is_right_cross(i)
                saccade_squarenum(iSac) = 6;
            else
                saccade_squarenum(iSac) = 0;
            end
            if any(saccade_squarenum(fixation_saccades(1:j-1))==saccade_squarenum(iSac)) % saccade within fixcross
                saccade_class(iSac) = Constants.WITHIN;
            else % saccade to fixcross
                saccade_class(iSac) = Constants.FIXCROSS;
            end            
        else % during the trial, but not to an object
            saccade_class(iSac) = Constants.OTHER;                                  
        end            
    end
    
    % Check for saccades to squares
    trial_saccades = find(times>sq.trial.start_time(i) & times<=sq.trial.end_time(i));
    saccade_trialnum(trial_saccades) = i;
    for j=1:numel(trial_saccades)
        iSac = trial_saccades(j);
        for k=1:numel(Constants.SQUARE_X)
            dx = Constants.SQUARE_X(k)-positions(iSac,1);
            dy = Constants.SQUARE_Y(k)-positions(iSac,2);
            disto = sqrt(dx*dx + dy*dy);
            if disto<sq.pixel_threshold % if saccade is to this square
                saccade_distance(iSac) = disto;
                saccade_squarenum(iSac) = k;
                if (sq.trial.is_right_cross(i) && any(saccade_squarenum(trial_saccades(1:j-1))<saccade_squarenum(iSac))) ...
                        || (~sq.trial.is_right_cross(i) && any(saccade_squarenum(trial_saccades(1:j-1))>saccade_squarenum(iSac))) % if we've already been to the next square
                    saccade_class(iSac) = Constants.BACKWARD;
                elseif any(saccade_squarenum(trial_saccades(1:j-1))==saccade_squarenum(iSac)) % if we've been to this square already 
                    saccade_class(iSac) = Constants.WITHIN;
                elseif ~sq.trial.is_target_color(i,k)
                    saccade_class(iSac) = Constants.DISTRACTOR;
                elseif sq.trial.target_squares_sofar(i,k)==1
                    saccade_class(iSac) = Constants.INTEGRATION;
                elseif sq.trial.target_squares_sofar(i,k)==2
                    saccade_class(iSac) = Constants.COMPLETION;
                else
                    saccade_class(iSac) = Constants.EXTRA;                    
                end
            end
        end
        
        % Check for saccades to end circle
        if isnan(saccade_squarenum(iSac))
            if ~sq.trial.is_right_cross(i)
                dx = Constants.RIGHTCROSS_X-positions(iSac,1);
                dy = Constants.RIGHTCROSS_Y-positions(iSac,2);
            else
                dx = Constants.LEFTCROSS_X-positions(iSac,1);
                dy = Constants.LEFTCROSS_Y-positions(iSac,2);
            end
            disto = sqrt(dx*dx + dy*dy);
            if disto<sq.pixel_threshold
                saccade_distance(iSac) = disto;
                if sq.trial.is_right_cross(i)
                    saccade_squarenum(iSac) = 0;
                else
                    saccade_squarenum(iSac) = 6;
                end
                if any(saccade_squarenum(trial_saccades(1:j-1))==saccade_squarenum(iSac)) % saccade within endcircle
                    saccade_class(iSac) = Constants.WITHIN;
                else % saccade to endcircle
                    saccade_class(iSac) = Constants.ENDCIRCLE;
                end
            else % during the trial, but not to an object
                saccade_class(iSac) = Constants.OTHER;                                  
            end
        end
    end
    
    % Check for saccades to end circle
    circle_saccades = find(times>sq.trial.end_time(i) & times<=sq.trial.circle_time(i));
    saccade_trialnum(circle_saccades) = i;
    for j=1:numel(circle_saccades)
        iSac = circle_saccades(j);
        if sq.trial.is_right_cross(i)
            dx = Constants.LEFTCROSS_X-positions(iSac,1);
            dy = Constants.LEFTCROSS_Y-positions(iSac,2);
        else
            dx = Constants.RIGHTCROSS_X-positions(iSac,1);
            dy = Constants.RIGHTCROSS_Y-positions(iSac,2);
        end
        disto = sqrt(dx*dx + dy*dy);
        if disto<sq.pixel_threshold
            saccade_distance(iSac) = disto;
            if sq.trial.is_right_cross(i)
                saccade_squarenum(iSac) = 0;
            else
                saccade_squarenum(iSac) = 6;
            end
            if any(saccade_squarenum([trial_saccades; circle_saccades(1:j-1)])==saccade_squarenum(iSac)) % saccade within endcircle
                saccade_class(iSac) = Constants.WITHIN;
            else % saccade to fixcross
                saccade_class(iSac) = Constants.ENDCIRCLE;
            end            
        else % during the trial, but not to an object
            saccade_class(iSac) = Constants.OTHER;                                  
        end            
    end

    
end