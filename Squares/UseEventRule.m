function [elTimes, codes, elWeights, iEvents] = UseEventRule(y,rule)

% [elTimes, codes, elWeights, iEvents] = UseEventRule(y,rule)
% 
% INPUTS:
% -y is a vector of squares data structs (loaded using LoadBehaviorData)
% -rule is a string indicating which type of event you want to find
%
% OUTPUTS:
% -elTimes is an n-element cell array in which elTimes{i} contains the
% eyelink times of the requested events in session i.
% -codes is an n-element cell array in which codes{i} contains a cell array 
% of codes the same length as elTimes{i}.  Currently, this means just an
% array of the rule over and over.
% -elWeights is an n-element cell array in which elWeights{i} contains the
% weights that each event in elTimes{i} should be given in a GLM regressor.
% -iEvents is the indices of the saccades or events in the data struct that 
% match this event rule.
%
% Created 11/28/11 by DJ.
% Updated 12/28/11 by DJ - added options Extra, Irrel, Circle-T, Circle-D,
%    Cross, Errant.
% Updated 1/10/12 by DJ - added combination events for GLM contrasts
% Updated 2/8/12 by DJ - added more SacStart options, changed a few names
% Updated 3/22/12 by DJ - added SqNumX options, where X is 1-5.
% Updated 6/15/12 by DJ - added LogOddsTarg/Dist/etc. options
% Updated 6/25/12 by DJ - added elWeights output. TO DO: TEST!!!
% Updated 7/16/12 by DJ - FIXED SqNum regressors (to disallow repeated 
%    saccades), added iEvents output.
% Updated 7/18/12 by DJ - FIXED SqNum regressor bug (reversing Sq1-Sq5!)
% Updated 7/20/12 by DJ - worked on DiffLogOdds suite of rules/weights
% Updated 7/23/12 by DJ - completed DiffLogOdds suite of rules/weights
% Updated 8/10/12 by DJ - added LSqNum and RSqNum (bit of a hack!)
% Updated 8/13/12 by DJ - added Circle, L/RCircle, L/RSaccade
% Updated 1/23/13 by DJ - added New-Dist-0T, New-Dist-1T, Incompl,
%   TrialStart-L/R, TrialEnd-L/R
% Updated 1/24/13 by DJ - added New-Integ and New-Irrel, fixed 
% (Diff)LogOddsIrrel to include SqNum5 distractors
% Updated 1/25/13 by DJ - added TrialStart/TrialEnd rules
% Updated 1/28/13 by DJ - fixed TrialStart/TrialEnd bug
% Updated 2/13/13 by DJ - added v2pt2 rules (debugged 2/19/13)
% Updated 5/17/13 by DJ - added sf3 compatibility
% Updated 3/10/14 by DJ - added "new" labels for sf, sf3
% Updated 3/11/14 by DJ - added pT_{0/2}-style labels, column vec check
% Updated 7/30/14 by DJ - allow rules starting with pT_{0/2}-style labels
%  and including ramps 
% Updated 1/28/15 by DJ - allow rules like p_D{x/3}

% Set up
nSessions = numel(y);
Constants = GetSquaresConstants;
elTimes = cell(1,nSessions);
codes = cell(1,nSessions);
elWeights = cell(1,nSessions);
iEvents = cell(1,nSessions);
% load logOdds variable if needed
if ~isempty(strfind(rule,'LogOdds'))
    logOddsRatio = GetLogOddsRatio(y);
end

% Main loop
for i=1:nSessions
    
    % analysis 3pt0 options
    if ~isempty(strfind(rule,'/3}')) % SquaresFix3
        [type,sqNum] = GetSquareTypes(y(i),'sf3',1); 
        if ~isempty(strfind(rule,'{x/'))
            iEvents{i} = find(strncmp(rule,type,2)); % pD or pT
        else
            iEvents{i} = find(strncmp(rule,type,find(rule=='}')));
        end
        elTimes{i} = y(i).trial.square_time(iEvents{i});
        
    elseif ~isempty(strfind(rule,'/2}')) && rule(1)=='p'  % SquaresFix
        [type,sqNum] = GetSquareTypes(y(i),'sf',1);        
        if ~isempty(strfind(rule,'{x/'))
            iEvents{i} = find(strncmp(rule,type,2)); % pD or pT
        else
            iEvents{i} = find(strncmp(rule,type,find(rule=='}')));
        end
        elTimes{i} = y(i).trial.square_time(iEvents{i});
    
    elseif ~isempty(strfind(rule,'/2}')) && rule(1)=='a'  % Squares
        [type,sqNum] = GetSquareTypes(y(i),'sq',1); 
        % get index of 1st saccade to each relevant square
        if ~isempty(strfind(rule,'{x/'))
            [rows, cols] = find(strncmp(rule,type,2)); % aD or aT
        else
            [rows, cols] = find(strncmp(rule,type,find(rule=='}')));
        end
        iEvents{i} = nan(1,length(rows));
        for j=1:length(rows) 
            % get first saccade to this trialnum, squarenum combo
            iFirst2sq = find(y(i).saccade.trialnum==rows(j) & y(i).saccade.squarenum==cols(j), 1);
            if ~isempty(iFirst2sq)
                iEvents{i}(j) = iFirst2sq;
            end
        end
        % remove the squares with no saccades to them
        rows = rows(~isnan(iEvents{i}));
        cols = cols(~isnan(iEvents{i}));
        iEvents{i} = iEvents{i}(~isnan(iEvents{i}));
        elTimes{i} = y(i).saccade.end_time(iEvents{i});
        
    
    % SquaresFix v1pt0 options (this is v1.0 of analysis, not experiment)
    elseif ~isempty(strfind(rule,'sf3-'));
        if isempty(strfind(rule,'new-'));   % e.g. sf3-
            useNew = false;
            newrule = rule(5:end);
        else    % e.g. new-sf3-
            useNew = true;
            newrule = rule(8:end); 
        end
        [type, sqNum] = GetSquareTypes(y(i),'sf3',useNew);
        nTrials = numel(y(i).trial.start_time);
        if strcmp(newrule,'FixCross')
            iEvents{i} = 1:nTrials;
            elTimes{i} = y(i).trial.fix_time(iEvents{i}); 
        elseif strcmp(newrule,'Circle')
            iEvents{i} = 1:nTrials;
            elTimes{i} = y(i).trial.end_time(iEvents{i}); 
        elseif strcmp(newrule,'Square')
            iEvents{i} = 1:numel(type);
            elTimes{i} = y(i).trial.square_time(iEvents{i}); 
        elseif strncmp(newrule,'SqNum',5)            
            num = str2double(newrule(6)); % 1-5
            iEvents{i} = find(sqNum==num);
            elTimes{i} = y(i).trial.square_time(iEvents{i}); 
        else            
            iEvents{i} = find(strcmp(newrule,type));
            elTimes{i} = y(i).trial.square_time(iEvents{i}); 
        end 
    
    elseif ~isempty(strfind(rule,'sf-'));
        if isempty(strfind(rule,'new-'));   % e.g. sf-
            useNew = false;
            newrule = rule(4:end);
        else    % e.g. new-sf-
            useNew = true;
            newrule = rule(7:end); 
        end
        [type, sqNum] = GetSquareTypes(y(i),'sf',useNew);
        nTrials = numel(y(i).trial.start_time);
        if strcmp(newrule,'FixCross')
            iEvents{i} = 1:nTrials;
            elTimes{i} = y(i).trial.fix_time(iEvents{i}); 
        elseif strcmp(newrule,'Circle')
            iEvents{i} = 1:nTrials;
            elTimes{i} = y(i).trial.end_time(iEvents{i}); 
        elseif strcmp(newrule,'Square')
            iEvents{i} = 1:numel(type);
            elTimes{i} = y(i).trial.square_time(iEvents{i}); 
        elseif strncmp(newrule,'SqNum',5)            
            num = str2double(newrule(6)); % 1-5
            iEvents{i} = find(sqNum==num);
            elTimes{i} = y(i).trial.square_time(iEvents{i}); 
        else            
            iEvents{i} = find(strcmp(newrule,type));
            elTimes{i} = y(i).trial.square_time(iEvents{i}); 
        end    
        
    elseif ~isempty(strfind(rule,'sq-'));
        if isempty(strfind(rule,'new-'));   % e.g. sq-
            useNew = false;
            newrule = rule(4:end);
        else    % e.g. new-sq-
            useNew = true;
            newrule = rule(7:end); 
        end
        [type, sqNum] = GetSquareTypes(y(i),'sq',useNew);
        nTrials = numel(y(i).trial.start_time);
        if strcmp(newrule,'FixCross')
            iEvents{i} = 1:nTrials;
            elTimes{i} = y(i).trial.fix_time(iEvents{i}); 
        elseif strcmp(newrule,'Circle')
            iEvents{i} = 1:nTrials;
            elTimes{i} = y(i).trial.end_time(iEvents{i}); 
        elseif strcmp(newrule,'Square')
            iEvents{i} = 1:numel(type);
            elTimes{i} = y(i).trial.square_time(iEvents{i});
        elseif strncmp(newrule,'SqNum',5)            
            num = str2double(newrule(6)); % 1-5
            iEvents{i} = find(sqNum==num);
            elTimes{i} = y(i).trial.square_time(iEvents{i}); 
        else            
            iEvents{i} = find(strcmp(newrule,type));
            elTimes{i} = y(i).trial.square_time(iEvents{i}); 
        end

        
        
    % Blink events (for SquaresFix datasets)
    elseif strcmp(rule,'BlinkStart')
        iEvents{i} = 1:numel(y(i).blink.start_time);
        elTimes{i} = y(i).blink.start_time(iEvents{i}); 
    elseif strcmp(rule,'BlinkEnd')
        iEvents{i} = 1:numel(y(i).blink.end_time);
        elTimes{i} = y(i).blink.end_time(iEvents{i}); 
    
    % Squares Version 2.2 options
    elseif ~isempty(strfind(rule,'-v2pt2'))
        
        if strcmp(rule,'MiniSaccade-v2pt2') % corrective saccade            
            iSacToSq = find(y(i).saccade.squarenum>0 & y(i).saccade.squarenum<6);
            dist = sqrt(sum((y(i).saccade.end_position(iSacToSq,:)-y(i).saccade.start_position(iSacToSq,:)).^2, 2));
            iEvents{i} = iSacToSq(dist<=50); % 1/3 of ditance between squares                
        elseif strcmp(rule,'LSaccade-v2pt2')
            iSacToSqOrCircle = find(y(i).saccade.squarenum>=0 & y(i).saccade.squarenum<6); % include 0 for leftward trials' circle
            iLSacToSqOrCircle = iSacToSqOrCircle(y(i).trial.is_right_cross(y(i).saccade.trialnum(iSacToSqOrCircle))); % crop to leftward trials only
            dist = sqrt(sum((y(i).saccade.end_position(iLSacToSqOrCircle,:)-y(i).saccade.start_position(iLSacToSqOrCircle,:)).^2, 2));
            iEvents{i} = iLSacToSqOrCircle(dist>50); % 1/3 of ditance between squares
        elseif strcmp(rule,'RSaccade-v2pt2')
            iSacToSqOrCircle = find(y(i).saccade.squarenum>0 & y(i).saccade.squarenum<=6); % include 6 for rightward trials' circle
            iRSacToSqOrCircle = iSacToSqOrCircle(~y(i).trial.is_right_cross(y(i).saccade.trialnum(iSacToSqOrCircle))); % crop to rightward trials only
            dist = sqrt(sum((y(i).saccade.end_position(iRSacToSqOrCircle,:)-y(i).saccade.start_position(iRSacToSqOrCircle,:)).^2, 2));
            iEvents{i} = iRSacToSqOrCircle(dist>50); % 1/3 of ditance between squares
            
        elseif strcmp(rule(1:5),'SqNum') % SqNum rules
            % Get longest-lingering saccade on each square
            [iSaccades,~,sqNum] = FindSaccadesToSquares(y(i),'longest'); % make 2nd argument 0 for 1st saccade, inf for last saccade
            num = str2double(rule(6)); % 1-5
            iEvents{i} = iSaccades(sqNum==num);
        else % Type rules
            % Get longest-lingering saccade on each square
            [iSaccades,type] = FindSaccadesToSquares(y(i),'longest'); % make 2nd argument 0 for 1st saccade, inf for last saccade
            newrule = rule(1:end-6);       
            iEvents{i} = iSaccades(strcmp(newrule,type));
        end
        % Remove NaN's
        iEvents{i} = iEvents{i}(~isnan(iEvents{i}));
        % fill in times
        elTimes{i} = y(i).saccade.end_time(iEvents{i});
        
    % Version 2.1 and before   
    else
    
        switch rule   
            
            

            case 'SacStartInteg'
                iEvents{i} = find(y(i).saccade.class==Constants.INTEGRATION);
                elTimes{i} = y(i).saccade.start_time(iEvents{i});
            case {'SacEndInteg' 'Integ'}
                iEvents{i} = find(y(i).saccade.class==Constants.INTEGRATION);
                elTimes{i} = y(i).saccade.end_time(iEvents{i});   
            case {'New-Integ','LogOddsInteg' 'DiffLogOddsInteg'} % exclude saccades where there's no chance of a complete trial (square 5)
                SqNum = y(i).saccade.squarenum;
                rightys = find(y(i).trial.is_right_cross);
                SqNum(ismember(y(i).saccade.trialnum,rightys)) = 6-SqNum(ismember(y(i).saccade.trialnum,rightys));
                iEvents{i} = find(y(i).saccade.class==Constants.INTEGRATION & SqNum<5);
                elTimes{i} = y(i).saccade.end_time(iEvents{i});
            case 'SacStartCompl'
                iEvents{i} = find(y(i).saccade.class==Constants.COMPLETION);
                elTimes{i} = y(i).saccade.start_time(iEvents{i});
            case {'SacEndCompl' 'Compl' 'LogOddsCompl' 'DiffLogOddsCompl'}
                iEvents{i} = find(y(i).saccade.class==Constants.COMPLETION);
                elTimes{i} = y(i).saccade.end_time(iEvents{i});                       
            case 'SacStartExtra'
                iEvents{i} = find(y(i).saccade.class==Constants.EXTRA);
                elTimes{i} = y(i).saccade.start_time(iEvents{i});
            case {'SacEndExtra' 'Extra'}
                iEvents{i} = find(y(i).saccade.class==Constants.EXTRA);
                elTimes{i} = y(i).saccade.end_time(iEvents{i});
            case 'SacStartErrant' % saccade to something other than the next (or current) square
                iEvents{i} = find(y(i).saccade.class==Constants.OTHER | y(i).saccade.class==Constants.BACKWARD);
                elTimes{i} = y(i).saccade.start_time(iEvents{i});
            case {'SacEndErrant' 'Errant'}
                iEvents{i} = find(y(i).saccade.class==Constants.OTHER | y(i).saccade.class==Constants.BACKWARD);
                elTimes{i} = y(i).saccade.end_time(iEvents{i});    
            case 'SacStartCross' % saccade to the fixation cross
                iEvents{i} = find(y(i).saccade.class==Constants.FIXCROSS);
                elTimes{i} = y(i).saccade.start_time(iEvents{i});
            case {'SacEndCross' 'Cross'}
                iEvents{i} = find(y(i).saccade.class==Constants.FIXCROSS);
                elTimes{i} = y(i).saccade.end_time(iEvents{i});
            case 'SacStartCircle-D' % saccade to the end circle on a distractor trial
                okSaccades = find(y(i).saccade.class==Constants.ENDCIRCLE & ~isnan(y(i).saccade.trialnum));
                isTargetTrial = y(i).trial.is_target_trial(y(i).saccade.trialnum(okSaccades));
                iEvents{i} = okSaccades(~isTargetTrial);
                elTimes{i} = y(i).saccade.start_time(iEvents{i});
            case {'SacEndCircle-D' 'Circle-D'}
                okSaccades = find(y(i).saccade.class==Constants.ENDCIRCLE & ~isnan(y(i).saccade.trialnum));
                isTargetTrial = y(i).trial.is_target_trial(y(i).saccade.trialnum(okSaccades));
                iEvents{i} = okSaccades(~isTargetTrial);
                elTimes{i} = y(i).saccade.end_time(iEvents{i});
            case 'SacStartCircle-T' % saccade to the end circle on a target trial
                okSaccades = find(y(i).saccade.class==Constants.ENDCIRCLE & ~isnan(y(i).saccade.trialnum));
                isTargetTrial = y(i).trial.is_target_trial(y(i).saccade.trialnum(okSaccades));
                iEvents{i} = okSaccades(isTargetTrial);
                elTimes{i} = y(i).saccade.start_time(iEvents{i});
            case {'SacEndCircle-T' 'Circle-T'}
                okSaccades = find(y(i).saccade.class==Constants.ENDCIRCLE & ~isnan(y(i).saccade.trialnum));
                isTargetTrial = y(i).trial.is_target_trial(y(i).saccade.trialnum(okSaccades));
                iEvents{i} = okSaccades(isTargetTrial);
                elTimes{i} = y(i).saccade.end_time(iEvents{i});
            case 'SacStartCircle'
                iEvents{i} = find(y(i).saccade.class==Constants.ENDCIRCLE & ~isnan(y(i).saccade.trialnum));            
                elTimes{i} = y(i).saccade.start_time(iEvents{i});
            case {'SacEndCircle' 'Circle'}
                iEvents{i} = find(y(i).saccade.class==Constants.ENDCIRCLE & ~isnan(y(i).saccade.trialnum));            
                elTimes{i} = y(i).saccade.end_time(iEvents{i});
            case 'LCircle' % leftward saccade to a circle
                okSaccades = find(y(i).saccade.class==Constants.ENDCIRCLE & ~isnan(y(i).saccade.trialnum));            
                isLeftTrial = y(i).trial.is_right_cross(y(i).saccade.trialnum(okSaccades));
                iEvents{i} = okSaccades(isLeftTrial);
                elTimes{i} = y(i).saccade.end_time(iEvents{i});
            case 'RCircle' % rightward saccade to a circle
                okSaccades = find(y(i).saccade.class==Constants.ENDCIRCLE & ~isnan(y(i).saccade.trialnum));            
                isLeftTrial = ~y(i).trial.is_right_cross(y(i).saccade.trialnum(okSaccades));
                iEvents{i} = okSaccades(isLeftTrial);
                elTimes{i} = y(i).saccade.end_time(iEvents{i});
            case {'SacStartIrrel' 'SacEndIrrel' 'Irrel' 'New-Irrel' 'LogOddsIrrel' 'DiffLogOddsIrrel'} % Irrel = IRRELEVANT. saccade to any square after the trial is over, whether or not it's a target
                % get number of targets seen so far - include Integration saccades for LogOddsIrrel
                okSaccades = find((y(i).saccade.class==Constants.DISTRACTOR | y(i).saccade.class==Constants.EXTRA | y(i).saccade.class==Constants.INTEGRATION) & y(i).saccade.squarenum>0 & y(i).saccade.squarenum<6);
                targets_sofar = nan(size(okSaccades));
                for j=1:numel(okSaccades)
                    targets_sofar(j) = y(i).trial.target_squares_sofar(y(i).saccade.trialnum(okSaccades(j)),y(i).saccade.squarenum(okSaccades(j)));
                end            
                switch rule
                    case 'SacStartIrrel'          
                        iEvents{i} = okSaccades(targets_sofar>=2);
                        elTimes{i} = y(i).saccade.start_time(iEvents{i});
                    case {'SacEndIrrel' 'Irrel'}
                        iEvents{i} = okSaccades(targets_sofar>=2);
                        elTimes{i} = y(i).saccade.end_time(iEvents{i});
                    case {'New-Irrel','LogOddsIrrel' 'DiffLogOddsIrrel'}
                        SqNum = y(i).saccade.squarenum;
                        rightys = find(y(i).trial.is_right_cross);
                        SqNum(ismember(y(i).saccade.trialnum,rightys)) = 6-SqNum(ismember(y(i).saccade.trialnum,rightys));
    %                     iEvents{i} = okSaccades(targets_sofar>=2 | (y(i).saccade.class(okSaccades)==Constants.INTEGRATION & SqNum(okSaccades)==5));
                        iEvents{i} = okSaccades(targets_sofar>=2 | (targets_sofar==0 & SqNum(okSaccades)==5));
                        elTimes{i} = y(i).saccade.end_time(iEvents{i});
                end             

            % COMBINATION EVENTS FOR GLM CONTRASTS
            case {'Preparation' 'Prep'}
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'Dist-0T');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'Integ');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});            
            case {'Anticipation' 'Antic'}
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'Dist-1T');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'Compl');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case {'Distractor' 'Dist'}
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'Dist-0T');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'Dist-1T');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case 'LogOddsDist'
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'LogOddsDist-0T');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'LogOddsDist-1T');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case {'Target' 'Targ'}
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'Integ');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'Compl');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case 'LogOddsTarg'
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'LogOddsInteg');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'LogOddsCompl');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case 'SacStartSquare'
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'SacStartDist-0T');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'SacStartInteg');
                [tempevents(3),~,~,temp_iEvents(3)] = UseEventRule(y(i),'SacStartDist-1T');            
                [tempevents(4),~,~,temp_iEvents(4)] = UseEventRule(y(i),'SacStartCompl');
                [tempevents(5),~,~,temp_iEvents(5)] = UseEventRule(y(i),'SacStartIrrel');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case {'Square'}
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'Dist-0T');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'Integ');
                [tempevents(3),~,~,temp_iEvents(3)] = UseEventRule(y(i),'Dist-1T');            
                [tempevents(4),~,~,temp_iEvents(4)] = UseEventRule(y(i),'Compl');
                [tempevents(5),~,~,temp_iEvents(5)] = UseEventRule(y(i),'Irrel');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case {'LSquare'}
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'LSqNum1');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'LSqNum2');
                [tempevents(3),~,~,temp_iEvents(3)] = UseEventRule(y(i),'LSqNum3');            
                [tempevents(4),~,~,temp_iEvents(4)] = UseEventRule(y(i),'LSqNum4');
                [tempevents(5),~,~,temp_iEvents(5)] = UseEventRule(y(i),'LSqNum5');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case {'RSquare'}
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'RSqNum1');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'RSqNum2');
                [tempevents(3),~,~,temp_iEvents(3)] = UseEventRule(y(i),'RSqNum3');            
                [tempevents(4),~,~,temp_iEvents(4)] = UseEventRule(y(i),'RSqNum4');
                [tempevents(5),~,~,temp_iEvents(5)] = UseEventRule(y(i),'RSqNum5');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case 'Sq1-4'
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'SqNum1');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'SqNum2');
                [tempevents(3),~,~,temp_iEvents(3)] = UseEventRule(y(i),'SqNum3');            
                [tempevents(4),~,~,temp_iEvents(4)] = UseEventRule(y(i),'SqNum4');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
%             case 'SacStartCircle'
%                 [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'SacStartCircle-D');
%                 [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'SacStartCircle-T');
%                 iEvents{i} = cat(1,temp_iEvents{:});
%                 elTimes{i} = cat(1,tempevents{:});
%             case {'SacEndCircle', 'Circle'}
%                 [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'Circle-D');
%                 [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'Circle-T');
%                 iEvents{i} = cat(1,temp_iEvents{:});
%                 elTimes{i} = cat(1,tempevents{:});
            case 'StartSaccade'
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'SacStartSquare');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'SacStartCircle');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case 'Saccade'
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'Square');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'Circle');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case 'LSaccade'
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'LSquare');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'LCircle');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case 'RSaccade'
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'RSquare');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'RCircle');
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});
            case 'Button'
                iEvents{i} = find(~isnan(y(i).trial.response_time)); % exclude trials with no response
                elTimes{i} = y(i).trial.response_time(iEvents{i}); 
            case 'DiffLogOddsSquare'
                [tempevents(1),~,~,temp_iEvents(1)] = UseEventRule(y(i),'DiffLogOddsDist-0T');
                [tempevents(2),~,~,temp_iEvents(2)] = UseEventRule(y(i),'DiffLogOddsDist-1T');
                [tempevents(3),~,~,temp_iEvents(3)] = UseEventRule(y(i),'DiffLogOddsInteg');            
                iEvents{i} = cat(1,temp_iEvents{:});
                elTimes{i} = cat(1,tempevents{:});

            % TRIAL START/END EVENTS
            case 'TrialStart'
                iEvents{i} = 1:numel(y(i).trial.start_time);
                elTimes{i} = y(i).trial.start_time(iEvents{i});
            case 'TrialStart-L'
                iEvents{i} = find(~y(i).trial.is_right_cross);
                elTimes{i} = y(i).trial.start_time(iEvents{i});
            case 'TrialStart-R'
                iEvents{i} = find(y(i).trial.is_right_cross);
                elTimes{i} = y(i).trial.start_time(iEvents{i});
            case 'TrialEnd'
                iEvents{i} = 1:numel(y(i).trial.start_time);
                elTimes{i} = y(i).trial.end_time(iEvents{i});
            case 'TrialEnd-L'
                iEvents{i} = find(~y(i).trial.is_right_cross);
                elTimes{i} = y(i).trial.end_time(iEvents{i});
            case 'TrialEnd-R'
                iEvents{i} = find(y(i).trial.is_right_cross);
                elTimes{i} = y(i).trial.end_time(iEvents{i});

            otherwise

                if strfind(rule,'RSqNum') % Square number, leftward saccade
    %                 nSquares = size(y(i).trial.is_target_color,2);
                    sqnum = str2double(rule(end));               
                    okSaccades = find(~isnan(y(i).saccade.trialnum) & ismember(y(i).saccade.class,...
                        [Constants.DISTRACTOR, Constants.INTEGRATION, Constants.COMPLETION, Constants.EXTRA]));
                    iEvents{i} = okSaccades( (y(i).saccade.squarenum(okSaccades)==sqnum & y(i).trial.is_right_cross(y(i).saccade.trialnum(okSaccades))==0) );
                    if strfind(rule,'Start') % e.g., SacStartSqNum1
                        elTimes{i} = y(i).saccade.start_time(iEvents{i});
                    else % e.g., SacEndSqNum2, SqNum3
                        elTimes{i} = y(i).saccade.end_time(iEvents{i});
                    end

                elseif strfind(rule,'LSqNum') % Square number
    %                 nSquares = size(y(i).trial.is_target_color,2);
                    sqnum = str2double(rule(end));               
                    okSaccades = find(~isnan(y(i).saccade.trialnum) & ismember(y(i).saccade.class,...
                        [Constants.DISTRACTOR, Constants.INTEGRATION, Constants.COMPLETION, Constants.EXTRA]));
                    iEvents{i} = okSaccades( (y(i).saccade.squarenum(okSaccades)==6-sqnum & y(i).trial.is_right_cross(y(i).saccade.trialnum(okSaccades))==1) );
                    if strfind(rule,'Start') % e.g., SacStartSqNum1
                        elTimes{i} = y(i).saccade.start_time(iEvents{i});
                    else % e.g., SacEndSqNum2, SqNum3
                        elTimes{i} = y(i).saccade.end_time(iEvents{i});
                    end

                elseif strfind(rule,'SqNum') % Square number
    %                 nSquares = size(y(i).trial.is_target_color,2);
                    sqnum = str2double(rule(end));               
                    okSaccades = find(~isnan(y(i).saccade.trialnum) & ismember(y(i).saccade.class,...
                        [Constants.DISTRACTOR, Constants.INTEGRATION, Constants.COMPLETION, Constants.EXTRA]));
                    iEvents{i} = okSaccades( (y(i).saccade.squarenum(okSaccades)==sqnum & y(i).trial.is_right_cross(y(i).saccade.trialnum(okSaccades))==0) ...
                        | (y(i).saccade.squarenum(okSaccades)==6-sqnum & y(i).trial.is_right_cross(y(i).saccade.trialnum(okSaccades))==1) );
                    if strfind(rule,'Start') % e.g., SacStartSqNum1
                        elTimes{i} = y(i).saccade.start_time(iEvents{i});
                    else % e.g., SacEndSqNum2, SqNum3
                        elTimes{i} = y(i).saccade.end_time(iEvents{i});
                    end

                else % Distractor events sorted by number of targets seen so far
                    % get number of targets seen so far
                    okSaccades = find(y(i).saccade.class==Constants.DISTRACTOR & y(i).saccade.squarenum>0 & y(i).saccade.squarenum<6);
                    targets_sofar = nan(size(okSaccades));
                    for j=1:numel(okSaccades)
                        targets_sofar(j) = y(i).trial.target_squares_sofar(y(i).saccade.trialnum(okSaccades(j)),y(i).saccade.squarenum(okSaccades(j)));
                    end            
                    switch rule
                        case 'SacStartDist-0T' 
                            iEvents{i} = okSaccades(targets_sofar==0);
                            elTimes{i} = y(i).saccade.start_time(iEvents{i});
                        case {'SacEndDist-0T' 'Dist-0T'}
                            iEvents{i} = okSaccades(targets_sofar==0);
                            elTimes{i} = y(i).saccade.end_time(iEvents{i});
                        case 'LogOddsDist-0T' % exclude saccades where there's no chance of a complete trial (square 5)
                            SqNum = y(i).saccade.squarenum;
                            rightys = find(y(i).trial.is_right_cross);
                            SqNum(ismember(y(i).saccade.trialnum,rightys)) = 6-SqNum(ismember(y(i).saccade.trialnum,rightys));
                            iEvents{i} = okSaccades(targets_sofar==0 & SqNum(okSaccades)<5);
                            elTimes{i} = y(i).saccade.end_time(iEvents{i});
                        case {'New-Dist-0T','DiffLogOddsDist-0T'} % exclude saccades where there's no chance of a complete trial after it (square 4-5)
                            SqNum = y(i).saccade.squarenum;
                            rightys = find(y(i).trial.is_right_cross);
                            SqNum(ismember(y(i).saccade.trialnum,rightys)) = 6-SqNum(ismember(y(i).saccade.trialnum,rightys));
                            iEvents{i} = okSaccades(targets_sofar==0 & SqNum(okSaccades)<4);
                            elTimes{i} = y(i).saccade.end_time(iEvents{i});
                        case 'SacStartDist-1T'
                            iEvents{i} = okSaccades(targets_sofar==1);
                            elTimes{i} = y(i).saccade.start_time(iEvents{i});
                        case {'SacEndDist-1T' 'Dist-1T' 'LogOddsDist-1T'}
                            iEvents{i} = okSaccades(targets_sofar==1);
                            elTimes{i} = y(i).saccade.end_time(iEvents{i});
                        case {'New-Dist-1T','DiffLogOddsDist-1T'} % exclude saccades where there's no chance of a complete trial after it (square 5)
                            SqNum = y(i).saccade.squarenum;
                            rightys = find(y(i).trial.is_right_cross);
                            SqNum(ismember(y(i).saccade.trialnum,rightys)) = 6-SqNum(ismember(y(i).saccade.trialnum,rightys));
                            iEvents{i} = okSaccades(targets_sofar==1 & SqNum(okSaccades)<5);
                            elTimes{i} = y(i).saccade.end_time(iEvents{i});                            
                        case 'SacStartDist-2T'
                            iEvents{i} = okSaccades(targets_sofar==2);
                            elTimes{i} = y(i).saccade.start_time(iEvents{i});
                        case {'SacEndDist-2T' 'Dist-2T'}
                            iEvents{i} = okSaccades(targets_sofar==2);
                            elTimes{i} = y(i).saccade.end_time(iEvents{i});
                        case {'Incompl','DiffLogOddsIncompl', 'DiffLogOddsIncompletion'} % trials that ensure an incomplete trial (D-1T on SqNum5, or D-0T on SqNum4)
                            SqNum = y(i).saccade.squarenum;
                            rightys = find(y(i).trial.is_right_cross);
                            SqNum(ismember(y(i).saccade.trialnum,rightys)) = 6-SqNum(ismember(y(i).saccade.trialnum,rightys));
                            iEvents{i} = okSaccades( (targets_sofar==1 & SqNum(okSaccades)==5) | (targets_sofar==0 & SqNum(okSaccades)==4));
                            elTimes{i} = y(i).saccade.end_time(iEvents{i});                            
                        otherwise
                            elTimes{i} = [];
                    end
                end
        end
    end
    
    
%     % Translate indices to event times
%     if ~isempty(strfind(rule,'Start'))
%         elTimes{i} = y(i).saccade.start_time(iEvents);
%     elseif ~isempty(strfind(rule,'Button'))
%         elTimes{i} = y(i).trial.response_time(iEvents);
%     else
%         elTimes{i} = y(i).saccade.end_time(iEvents);
%     end
            
    % Specify Weights
    switch rule
        case {'LogOddsInteg', 'LogOddsCompl', 'LogOddsDist-0T', 'LogOddsDist-1T'} % Non-irrelevant saccades with LogOdds in title
            elWeights{i} = nan(size(elTimes{i}));
            for j=1:numel(elWeights{i})
                elWeights{i}(j) = logOddsRatio{i}(y(i).saccade.trialnum(iEvents{i}(j)), y(i).saccade.squarenum(iEvents{i}(j)));
            end
        case {'DiffLogOddsInteg','DiffLogOddsDist-0T','DiffLogOddsDist-1T' 'DiffLogOddsSquare'} % Non-irrelevant-or-completion-or-incompletion saccades with DiffLogOdds in title        
            isRighty = y(i).trial.is_right_cross;
            diffLogOdds_0 = diff(logOddsRatio{i},1,2);
            diffLogOdds = nan(size(diffLogOdds_0,1),5);
            diffLogOdds(~isRighty,1:4) = diffLogOdds_0(~isRighty,:);            
            diffLogOdds(isRighty,2:5) = -diffLogOdds_0(isRighty,:); % R-L diff is the opposite of L-R diff
            % Find diffLogOdds of each saccade
            elWeights{i} = nan(size(elTimes{i}));
            for j=1:numel(elWeights{i})
                elWeights{i}(j) = diffLogOdds(y(i).saccade.trialnum(iEvents{i}(j)), y(i).saccade.squarenum(iEvents{i}(j)));
            end            
        otherwise
            bracket = strfind(rule,'}');
            if ~isempty(bracket) && length(rule)>bracket+1
                suffix = rule((bracket+2):end);
                switch suffix
                    case 'RampUp'
                        weights = [1 2 3 4 5];                        
                    case 'RampDown'
                        weights = [5 4 3 2 1];                        
                    case 'Peak'
                        weights = [1 2 3 2 1];                        
                    case 'Valley'
                        weights = [3 2 1 2 3];
                    otherwise
                        error('Suffix not recognized!')
                end                
                weights = weights/mean(weights); % normalize
                if rule(1)=='a' % active task
                    elWeights{i} = weights(sqNum(sub2ind(size(sqNum),rows,cols)))';
                else % passive tasks   
                    elWeights{i} = weights(sqNum(iEvents{i}))';
                end
            end
            if isempty(elWeights{i})
                elWeights{i} = ones(size(elTimes{i}));
            end
    end

    
    % Put events in chronological order (for combo event types)
    [elTimes{i}, order] = sort(elTimes{i},'ascend'); 
    elWeights{i} = elWeights{i}(order);
    iEvents{i} = iEvents{i}(order);
    % Define codes
    codes{i} = repmat({rule},size(elTimes{i}));
    
    if size(elTimes{i},2)>1 % FORCE COLUMN VECTORS
        elTimes{i} = elTimes{i}';
    end
    if size(codes{i},2)>1 % FORCE COLUMN VECTORS
        codes{i} = codes{i}';
    end
    if size(elWeights{i},2)>1 % FORCE COLUMN VECTORS
        elWeights{i} = elWeights{i}';
    end
end

