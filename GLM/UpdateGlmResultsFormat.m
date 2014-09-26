function results_new = UpdateGlmResultsFormat(results_old)

% Takes a set of results made from the old GLM scripts and puts it them the
% format of the new Gui.
%
% results_new = UpdateGlmResultsFormat(results_old)
%
% INPUT:
% -results_old is a vector of data structs loaded from a GLM.  If they're
% old, this program will update them.
%
% OUTPUT:
% -Rnew is the same as results_old, but with the format the same as
% what would come out of RunGlmGui.
%
% Created 7/20/12 by DJ.
% Updated 8/9/12 by DJ - added vThresh
% Updated 3/22/13 by DJ - added artifact_influence, 2-element influences.
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 3/19/14 by DJ - influence, artifact_influence, tResponse in cells
% Updated 8/14/14 by DJ - added lambda field
% Updated 8/20/14 by DJ - made lamda a cell
% Updated 8/27/14 by DJ - added demeanX and normalizeX

nResults = numel(results_old);

for i=1:nResults
    % make a single struct so we can change the fields without an error
    Rold = results_old(i);
    Rnew = results_old(i);

    if isfield(Rold,'NewEEG')
        Rnew.EEG = Rold.NewEEG;
        Rnew = rmfield(Rnew,'NewEEG');
    end

    if ~isfield(Rold,'iLevel');
        if ~isfield(Rold,'regressor_events')
            Rnew.regressor_events = {Rold.nuisance_events};
            Rnew.iLevel = 1;
        else
            Rnew.regressor_events = {Rold.nuisance_events, Rold.regressor_events};
            Rnew.iLevel = 2;
        end
    end
    
    if isfield(Rold,'nuisance_events');
        Rnew = rmfield(Rnew,'nuisance_events');
    end

    if ~isfield(Rold,'filenames');
        Rnew.filenames = {'Unknown'};
    end

    if ~isfield(Rold,'dataset')
        Rnew.dataset = Rnew.EEG.filename;
    end

    if ~isfield(Rold,'offset')
        Rnew.offset = NaN;
    end

    if ~isfield(Rold,'influence')
        Rnew.influence = 500;
    end

    if ~isfield(Rold,'stddev')
        Rnew.stddev = 0;
    end
    
    if ~isfield(Rold,'vThresh')
        Rnew.vThresh = Inf;
    end

    if ~isfield(Rold,'method')
        Rnew.method = 'mvregress';
    end

    if ~isfield(Rold,'trial_rej_rules')
        Rnew.trial_rej_rules = {'skipped_ends'  'skip'  'backward'};
    end

    if ~isfield(Rold,'artifact_events')
        Rnew.artifact_events = {'Button'  'BlinkEnd'  'BlinkStart'  'Errant'  'Cross'};
    end
    if ~isfield(Rold,'artifact_influence')
        Rnew.artifact_influence = Rnew.influence;
    end
    
    if ~iscell(Rnew.influence) && length(Rnew.influence)==1
        Rnew.influence = [-Rnew.influence, Rnew.influence];
    end
    if ~iscell(Rnew.artifact_influence) && length(Rnew.artifact_influence)==1
        Rnew.artifact_influence = [-Rnew.artifact_influence, Rnew.artifact_influence];
    end
    % new 4/30/13
    if ~iscell(Rold.responseFns) 
        Rnew.responseFns = cell(1,Rnew.iLevel);
        Rnew.responseFns{Rnew.iLevel} = Rold.responseFns;
    end
    % New 3/19/14
    if ~iscell(Rold.influence)
        Rnew.influence = repmat({Rold.influence},1,Rnew.iLevel);
    end
    if ~iscell(Rold.artifact_influence)
        Rnew.artifact_influence = repmat({Rold.artifact_influence},1,Rnew.iLevel);
    end
    if ~iscell(Rold.tResponse)
        Rnew.tResponse = repmat({Rold.tResponse},1,Rnew.iLevel);
    end
    % New 8/14/14
    if ~isfield(Rold,'lambda')
        Rnew.lambda = repmat({NaN},1,Rnew.iLevel);
    % New 8/20/14
    elseif ~iscell(Rold.lambda)
        Rnew.lambda = repmat({Rold.lambda},1,Rnew.iLevel);
    end
    % New 8/27/14
    if ~isfield(Rold,'normalizeX')
        Rnew.normalizeX = false;
    end
    if ~isfield(Rold,'demeanX')
        Rnew.demeanX = false;
    end
       
    % Add to output vector
    results_new(i) = Rnew;
end