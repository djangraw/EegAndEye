function weights = ApplySvdToMultiSubjects(master,results,tRange,responseNum,nComps,doPlot)

% ApplySvdToMultiSubjects(master,results,tRange,responseNum,nComps,doPlot)
%
% Created 1/5/12 by DJ.
% Updated 4/17/12 by DJ - added cell support, weights output
% Updated 7/18/12 by DJ - added nComps input
% Updated 7/20/12 by DJ - added compatibility for new RunGlmGui GLM's
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 10/2/13 by DJ - added doPlot input

if isempty(master)
    useMaster = 0;
else
    useMaster = 1;
end

if nargin<6 || isempty(doPlot)
    doPlot = 1;
end

% Update data structs
if useMaster % for backwards compatibility for before 7/19/12
    master = UpdateGlmResultsFormat(master);    
end
results = UpdateGlmResultsFormat(results);

if useMaster % Get weights from master, apply to results
    % Get list of all channels
    channels = {master.EEG.chanlocs.labels};
    for i=1:numel(results)
        channels = intersect(channels, {results(i).EEG.chanlocs.labels});
    end

    % Get indices of ok channels
    master.responseFns = master.responseFns{master.iLevel}(ismember({master.EEG.chanlocs.labels},channels),:,:);
    for i=1:numel(results)    
        results(i).responseFns = results(i).responseFns{results.iLevel}(ismember({results(i).EEG.chanlocs.labels},channels),:,:);
    end

    % Get and plot SVD results
    if doPlot
        figure(100); clf;
        MakeFigureTitle(master.EEG.filename(1:9),0);
        chanlocs = master.EEG.chanlocs(ismember({master.EEG.chanlocs.labels},channels));
    else
        chanlocs = [];
    end
    [weights, eigenvalues] = ApplySvdToGlmResults(master.responseFns(:,:,responseNum),master.tResponse,tRange,chanlocs,master.regressor_events{master.iLevel}(responseNum),nComps);
    
    if doPlot
        for i=1:numel(results)        
            figure(100+i); clf;
            MakeFigureTitle(results(i).EEG.filename(1:9),0);
            PlotSvdWeightsAndCourses(results(i).responseFns(:,:,responseNum),results(i).tResponse,weights,eigenvalues,chanlocs,results(i).regressor_events{results(i).iLevel}(responseNum),nComps);
        end
    end
else % Get weights from results, apply only to self
    weights = cell(1,numel(results));
    for i=1:numel(results)
        % extract results
        if isstruct(results)
            R = results(i);
        elseif iscell(results)
            R = results{i};
        end
        
        % Set up plot
        if doPlot
            figure(100+i); clf;
            MakeFigureTitle(R.EEG.filename(1:9),0);
            chanlocs = R.EEG.chanlocs;            
        else
            chanlocs = [];
        end
        
        % Get SVD weights
        weights{i} = ApplySvdToGlmResults(R.responseFns{R.iLevel}(:,:,responseNum),R.tResponse{R.iLevel},tRange,chanlocs,R.regressor_events{R.iLevel}(responseNum),nComps);
            
    end
end