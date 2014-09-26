function event_weights = GetEventWeights(y,rule)

% Returns the regressor weights for each of the requested events in a 
% squares experiment, for use in a GLM.
%
% event_weights = GetEventWeights(y,rule)
%
% INPUTS:
% -y is an n-element vector of behavioral data structs from the squares
% experiment.
% -rule is a string indicating the event type.
%
% OUTPUTS:
% -event_weights is a 1xn vector of cells, where event_weights{i} is a
% vector of event weights for the relevant events in session i.
%
% Created 6/18/12 by DJ based on UseEventRule.
%
% TO DO: Combine the two programs!

if isnumeric(y)
    y = loadBehaviorData(y);
end

nSessions = length(y);
event_weights = cell(1,nSessions);
logOddsRatio = GetLogOddsRatio(y);
Constants = GetSquaresConstants;

for i=1:nSessions
    % set up
    SqNum = y(i).saccade.squarenum;
    rightys = find(y(i).trial.is_right_cross);
    SqNum(ismember(y(i).saccade.trialnum,rightys)) = 6-SqNum(ismember(y(i).saccade.trialnum,rightys));
    % Get saccades
    switch rule
        case 'LogOddsTarg' % exclude saccades where there's no chance of a complete trial (no targets seen yet, square 5)            
            iSaccades = find(y(i).saccade.class==Constants.COMPLETION | (y(i).saccade.class==Constants.INTEGRATION & SqNum<5));        
        case 'LogOddsDist' % exclude saccades where there's no chance of a complete trial (no targets seen yet, square 5)
            % get number of targets seen so far
            okSaccades = find(y(i).saccade.class==Constants.DISTRACTOR & y(i).saccade.squarenum>0 & y(i).saccade.squarenum<6);
            targets_sofar = nan(size(okSaccades));
            for j=1:numel(okSaccades)
                targets_sofar(j) = y(i).trial.target_squares_sofar(y(i).saccade.trialnum(okSaccades(j)),y(i).saccade.squarenum(okSaccades(j)));
            end            
            iSaccades = (okSaccades(targets_sofar==1 | (targets_sofar==0 & SqNum(okSaccades)<5)));    
%         case 'LogOddsIrrel' % Irrel = IRRELEVANT. saccade to any square after the trial is over, whether or not it's a target
%             % get number of targets seen so far - include Integration saccades for LogOddsIrrel
%             okSaccades = find((y(i).saccade.class==Constants.DISTRACTOR | y(i).saccade.class==Constants.EXTRA | y(i).saccade.class==Constants.INTEGRATION) & y(i).saccade.squarenum>0 & y(i).saccade.squarenum<6);
%             targets_sofar = nan(size(okSaccades));
%             for j=1:numel(okSaccades)
%                 targets_sofar(j) = y(i).trial.target_squares_sofar(y(i).saccade.trialnum(okSaccades(j)),y(i).saccade.squarenum(okSaccades(j)));
%             end            
%             iSaccades = (okSaccades(targets_sofar>=2 | (targets_sofar==0 & y(i).saccade.squarenum(okSaccades)==5)));
        case 'DiffLogOddsSquare'
            % get number of targets seen so far
            okSaccades = find( (y(i).saccade.class==Constants.DISTRACTOR | y(i).saccade.class==Constants.INTEGRATION) & y(i).saccade.squarenum>0 & y(i).saccade.squarenum<6);
            targets_sofar = nan(size(okSaccades));
            for j=1:numel(okSaccades)
                targets_sofar(j) = y(i).trial.target_squares_sofar(y(i).saccade.trialnum(okSaccades(j)),y(i).saccade.squarenum(okSaccades(j)));
            end            
            iSaccades = (okSaccades(targets_sofar==1 | (targets_sofar==0 & SqNum(okSaccades)<5)));    
        case 'DiffLogOddsIncompletion'
            iSaccades = [];
        otherwise
            iSaccades = [];
    end
    elTimes{i} = y(i).saccade.end_time(iSaccades)';
    event_weights{i} = nan(size(elTimes{i}));
    for j=1:numel(event_weights{i})
        event_weights{i}(j) = logOddsRatio{i}(y(i).saccade.trialnum(iSaccades(j)), y(i).saccade.squarenum(iSaccades(j)));
    end
end