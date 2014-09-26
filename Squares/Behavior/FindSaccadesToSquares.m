function [iSaccades,types,sqNum,nones,multis] = FindSaccadesToSquares(sq,cutoff)

% This function chooses between multiple saccades to the same object based
% on a cutoff of fixation time (or the longest fixation time).
% 
% [iSaccades,types,sqNum] = FindSaccadesToSquares(sq,cutoff)
% [iSaccades,types,sqNum] = FindSaccadesToSquares(sq,'longest')
%
% INPUTS:
% - sq is a squares dataset imported by import_squares_data.m.
% - cutoff is the fixation time in ms below which a fixation doesn't count.
%  'The' saccade to each object is the first one above the cutoff. If no
%  saccades are above the cutoff, the last saccade to the square is used.
%
% OUTPUTS:
% - iSaccades is an Nx5 matrix (where N is the number of trials).
% iSaccades(i,j) is the index of the chosen saccade to square j in trial i.
% - types is an Nx5 cell array of strings, where types{i,j} is the square
% type of square j in trial i. These are 'v2pt2' events: the options are
% {'Dist-0T','Integ','Dist-1T','Compl','Incompl','Irrel'}.
% - SqNum is an Nx5 matrix where SqNum(i,j) is the number in the sequence 
% of square j in trial i.
%
% Created 2/12/13 by DJ.
% Updated 2/13/13 by DJ - SqNum input, plotting option

% Handle inputs
if nargin<2
    cutoff = 200;
end
if isequal(cutoff,'longest')
    fprintf('subject %d, session %d: Finding longest-fixation saccade to each square...\n',sq.subject,sq.session)
    uselongest = true;
else
    fprintf('subject %d, session %d: Finding first saccade over %g ms to each square...\n',sq.subject,sq.session, cutoff)    
    uselongest = false;
end
% Set up
nTrials = numel(sq.trial.start_time);
iSaccades = nan(nTrials,5);
types = cell(nTrials,5);
sqNum = zeros(nTrials,5);
% keep track of how many squares have 0 or >1 saccades over cutoff
nones = 0;
multis = 0;

for i=1:nTrials % trial number   
    for j=1:5 % square number (L-to-R)
        % Get SqNum (in sequence)
        if sq.trial.is_right_cross(i)
            sqNum(i,j) = 6-j; % square #5 is the first seen by the subject
        else
            sqNum(i,j) = j; % square #1 is the first seen by the subject
        end
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
        
            
        % Find saccades to this square
        iThisSq = find(sq.saccade.trialnum==i & sq.saccade.squarenum==j);
        if isempty(iThisSq)
            fprintf('Trial %d, Square %d: no fixations.\n',i,j);
        elseif numel(iThisSq)==1
            iSaccades(i,j) = iThisSq;
        else
            % Pick the winning sacccade!
            fixtime = sq.saccade.start_time(iThisSq+1) - sq.saccade.end_time(iThisSq);           
            ft{i,j} = fixtime;
            dist{i,j} = sqrt(sum((sq.saccade.end_position(iThisSq,:)-sq.saccade.start_position(iThisSq,:)).^2, 2));
            if uselongest % assume the important saccade started the longest fixation
                [~, iWinner] = max(fixtime);
                iSaccades(i,j) = iThisSq(iWinner);
            else % Find first saccade where the subject fixated for more than cutoff ms.            
                iWinner = find(fixtime>cutoff);
                if isempty(iWinner) % if no saccades exceeded threshold, use the last one to this square.
                    fprintf('Trial %d, Square %d: 0/%d fixations over %.0f ms... using last saccade.\n',i,j,numel(fixtime),cutoff);
                    nones = nones + 1;
                    iSaccades(i,j) = iThisSq(end);
                elseif numel(iWinner)==1 % if 1 saccade exceeded the cutoff, use it!
                    iSaccades(i,j) = iThisSq(iWinner);            
                else % if >1 saccade exceeded the cutoff, use the first one.
                    fprintf('Trial %d, Square %d: %d/%d fixations over %.0f ms... using first one.\n',i,j,numel(iWinner),numel(fixtime),cutoff);
                    multis = multis + 1;
                    iSaccades(i,j) = iThisSq(iWinner(1));
                end
            end
        end
    end
            
end

% Plot results
doplot = false;
if doplot
    iS = iSaccades(~isnan(iSaccades));
    f = sq.saccade.start_time(iS+1) - sq.saccade.end_time(iS);           
    d = sqrt(sum((sq.saccade.end_position(iS,:)-sq.saccade.start_position(iS,:)).^2, 2)); 

    cla; hold on;
    for i=1:numel(ft)
        plot(dist{i},ft{i},'.-')
        if ~isempty(dist{i})
            plot(dist{i}(end),ft{i}(end),'r.')
        end        
    end
    plot(d,f,'go'); % plot chosen saccades
    if ~ischar(cutoff) && ~isinf(cutoff)
        plot(get(gca,'xlim'),[cutoff cutoff],'k--','linewidth',2);
    end
    title(show_symbols(sq.eyeFilename));
    xlabel('distance (pixels');
    ylabel('fixation duration (ms)');
    axis([0 250 0 1000]);
    h = gca;
    MakeLegend({'b.-','r-','go'},{'saccades to same square','last saccade before leaving square','chosen saccade'});
    axes(h);
end
disp('done!')