function GetNumberMeaning(number)

% Displays the meaning of an event code in the Numbers struct.
% 
% GetNumberMeaning(number)
%
% INPUTS:
% - number is an integer used as an event code in EEGLAB.
% See GetNumbers.m for more details.
%
% Created 10/18/10 by DJ.

GetNumbers;

if number>Numbers.SACCADE_TO && number<Numbers.SACCADE_TO+50
    fprintf('SACCADE_TO object # %d\n', number-Numbers.SACCADE_TO);
elseif number>Numbers.ENTERS && number<Numbers.ENTERS+50
    fprintf('ENTER object # %d\n', number-Numbers.ENTERS); 
elseif number>Numbers.EXITS && number<Numbers.EXITS+50
    fprintf('EXIT object # %d\n', number-Numbers.EXITS); 
elseif number==Numbers.START_TRIAL+Numbers.STATIONARY
    fprintf('START Stationary Trial\n'); 
elseif number==Numbers.START_TRIAL+Numbers.MOVING
    fprintf('START Moving Trial\n'); 
elseif number==Numbers.START_TRIAL+Numbers.POPUP
    fprintf('START Popup Trial\n');     
else
    names = fieldnames(Numbers);
    for i=1:numel(names)
        if Numbers.(names{i})==number
            fprintf('%s event\n',names{i})
            return;
        end
    end
end
    