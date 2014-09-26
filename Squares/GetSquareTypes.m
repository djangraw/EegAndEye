function [types, sqNum] = GetSquareTypes(sq,prefix,useNew)

% Get each square's semantic type based on the task constraints.
%
% [types, sqNum] = GetSquareTypes(sq,prefix,useNew)
%
% INPUTS:
% - sq is a squares/sqfix/sf3 dataset struct.
% - prefix is a string indicating the dataset type (sq, sf, sf3).
% - useNew is a binary value indicating whether to use the new tagging
% system.
%
% OUTPUTS:
% - types is an NxM cell array of strings indicating the type of each 
% square in each trial ([N,M] is the size of sq.trial.is_target_color).
% - sqNum is an NxM double array indicating the number in the sequence of 
% each square number in each trial.
%
% Created 3/6/13 by DJ based on FindSaccadesToSquares.
% Updated 5/17/13 by DJ - added prefix input, sf3 compatibility
% Updated 5/20/13 by DJ - fixed Targ1 and Targ2 SqNum requirements
% Updated 3/10/14 by DJ - added useNew input (new labels, separate T+/D+)

if nargin<2 || isempty(prefix)
    prefix = 'sf';
end
if nargin<3 || isempty(useNew)
    useNew = false;
end

% Set up
nTrials = numel(sq.trial.start_time);
types = cell(nTrials,5);
sqNum = repmat(1:5,nTrials,1);
if isfield(sq.trial,'is_right_cross')
    sqNum(sq.trial.is_right_cross,:) = fliplr(sqNum(sq.trial.is_right_cross,:));
end

% Main loop
if useNew
    if strcmp(prefix,'sf3')
        for i=1:nTrials
            for j=1:5
                % Classify square        
                if sq.trial.is_target_color(i,j)
                    if sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)<4
                        types{i,j} = 'pT_{0/3}';
                    elseif sq.trial.target_squares_sofar(i,j)==2 && sqNum(i,j)<5
                        types{i,j} = 'pT_{1/3}';
                    elseif sq.trial.target_squares_sofar(i,j)==3
                        types{i,j} = 'pT^*_{2/3}';
                    elseif sq.trial.target_squares_sofar(i,j)>3
                        types{i,j} = 'pT_{+/3}'; % a target after a complete decision
                    else
                        types{i,j} = 'pT_{-/3}'; % a target after an incomplete decision
                    end
                else
                    if (sq.trial.target_squares_sofar(i,j)==0 && sqNum(i,j)==3) || (sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)==4) || (sq.trial.target_squares_sofar(i,j)==2 && sqNum(i,j)==5)
                        types{i,j} = 'pD^*_{-/3}';
                    elseif sq.trial.target_squares_sofar(i,j)==0 && sqNum(i,j)<3
                        types{i,j} = 'pD_{0/3}';
                    elseif sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)<4
                        types{i,j} = 'pD_{1/3}';
                    elseif sq.trial.target_squares_sofar(i,j)==2 && sqNum(i,j)<5
                        types{i,j} = 'pD_{2/3}';
                    elseif sq.trial.target_squares_sofar(i,j)>2
                        types{i,j} = 'pD_{+/3}';  % a distractor after a complete decision
                    else
                        types{i,j} = 'pD_{-/3}';  % a distractor after an incomplete decision
                    end
                end
            end
        end

    elseif strcmp(prefix,'sf')
        for i=1:nTrials
            for j=1:5
                % Classify square        
                if sq.trial.is_target_color(i,j)
                    if sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)~=5
                        types{i,j} = 'pT_{0/2}';
                    elseif sq.trial.target_squares_sofar(i,j)==2
                        types{i,j} = 'pT^*_{1/2}';
                    elseif sq.trial.target_squares_sofar(i,j)>2
                        types{i,j} = 'pT_{+/2}'; % a target after a complete decision                
                    else
                        types{i,j} = 'pT_{-/2}'; % a target after an incomplete decision
                    end
                else
                    if (sq.trial.target_squares_sofar(i,j)==0 && sqNum(i,j)==4) || (sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)==5)
                        types{i,j} = 'pD^*_{-/2}';
                    elseif sq.trial.target_squares_sofar(i,j)==0 && sqNum(i,j)<4
                        types{i,j} = 'pD_{0/2}';
                    elseif sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)<5
                        types{i,j} = 'pD_{1/2}';            
                    elseif sq.trial.target_squares_sofar(i,j)>1
                        types{i,j} = 'pD_{+/2}';  % a distractor after a complete decision
                    else
                        types{i,j} = 'pD_{-/2}';  % a distractor after an incomplete decision
                    end
                end
            end
        end

    elseif strcmp(prefix,'sq')
        for i=1:nTrials
            for j=1:5
                % Classify square        
                if sq.trial.is_target_color(i,j)
                    if sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)~=5
                        types{i,j} = 'aT_{0/2}';
                    elseif sq.trial.target_squares_sofar(i,j)==2
                        types{i,j} = 'aT^*_{1/2}';
                    elseif sq.trial.target_squares_sofar(i,j)>2
                        types{i,j} = 'aT_{+/2}'; % a target after a complete decision                
                    else
                        types{i,j} = 'aT_{-/2}'; % a target after an incomplete decision
                    end
                else
                    if (sq.trial.target_squares_sofar(i,j)==0 && sqNum(i,j)==4) || (sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)==5)
                        types{i,j} = 'aD^*_{-/2}';
                    elseif sq.trial.target_squares_sofar(i,j)==0 && sqNum(i,j)<4
                        types{i,j} = 'aD_{0/2}';
                    elseif sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)<5
                        types{i,j} = 'aD_{1/2}';            
                    elseif sq.trial.target_squares_sofar(i,j)>1
                        types{i,j} = 'aD_{+/2}';  % a distractor after a complete decision
                    else
                        types{i,j} = 'aD_{-/2}';  % a distractor after an incomplete decision
                    end
                end
            end
        end
        
    end
    
    
    
else %if ~useNew


    if strcmp(prefix,'sf3')
        for i=1:nTrials
            for j=1:5
                % Classify square        
                if sq.trial.is_target_color(i,j)
                    if sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)<4
                        types{i,j} = 'Targ1';
                    elseif sq.trial.target_squares_sofar(i,j)==2 && sqNum(i,j)<5
                        types{i,j} = 'Targ2';
                    elseif sq.trial.target_squares_sofar(i,j)==3
                        types{i,j} = 'Targ3';
                    else
                        types{i,j} = 'Irrel';
                    end
                else
                    if (sq.trial.target_squares_sofar(i,j)==0 && sqNum(i,j)==3) || (sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)==4) || (sq.trial.target_squares_sofar(i,j)==2 && sqNum(i,j)==5)
                        types{i,j} = 'Incompl';
                    elseif sq.trial.target_squares_sofar(i,j)==0 && sqNum(i,j)<3
                        types{i,j} = 'Dist-0T';
                    elseif sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)<4
                        types{i,j} = 'Dist-1T';
                    elseif sq.trial.target_squares_sofar(i,j)==2 && sqNum(i,j)<5
                        types{i,j} = 'Dist-2T';
                    else
                        types{i,j} = 'Irrel';  
                    end
                end
            end
        end

    elseif strcmp(prefix,'sf')
        for i=1:nTrials
            for j=1:5
                % Classify square        
                if sq.trial.is_target_color(i,j)
                    if sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)~=5
                        types{i,j} = 'Integ';
                    elseif sq.trial.target_squares_sofar(i,j)==2
                        types{i,j} = 'Compl';
                    else
                        types{i,j} = 'Irrel';                    
                    end
                else
                    if (sq.trial.target_squares_sofar(i,j)==0 && sqNum(i,j)==4) || (sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)==5)
                        types{i,j} = 'Incompl';
                    elseif sq.trial.target_squares_sofar(i,j)==0 && sqNum(i,j)<4
                        types{i,j} = 'Dist-0T';
                    elseif sq.trial.target_squares_sofar(i,j)==1 && sqNum(i,j)<5
                        types{i,j} = 'Dist-1T';            
                    else
                        types{i,j} = 'Irrel'; 
                    end
                end
            end
        end

    end
end