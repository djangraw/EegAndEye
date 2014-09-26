function [RF, meanRF] = PlotResponseFnsForThesis(R_cell,iEvents_cell,chansToPlot,iLevel_vec,addSquareRf)

% PlotResponseFnsForThesis(R_cell,iEvents_cell,chansToPlot,iLevel_vec,addSquareRf)
% 
% INPUTS:
% - R_cell is an n-element vector of cells containing vectors of structs of
% GLM results.
% - iEvents_cell is an n-element vector of cells containing the indices of
% the events you want to use for each experiment type.
% - chansToPlot is a cell array of strings indicating the channels whose
% response functions should be plotted. default: {'Fz';'Cz';'Pz';'Oz'}
% - iLevel_vec is an n-element vector of values representing the level of
% analysis you'd like to use (given a multi-level GLM). default: last.
% - addSquareRf is a binary value indicating whether you'd like to add the
% 'square' RF (the generic level 1 response) to each RF of interest. This
% will help visualize the changes in the context of the evoked potential.
%
% OUTPUTS:
% - RF is an n-element vector of cells, each containing the DxTxMxN matrix
% (where D=#chan, T=#samples, M=#events, and N=#subjects) representing the
% response function across all channels. 
% - meanRF is an n-element vecot of cells, each containing the DxTxM matrix
% doing the same, but averaged across subjects. 
%
% Created 3/13/14 by DJ.
% Updated 3/18/14 by DJ - added iLevel_cell input.
% Updated 3/19/14 by DJ - added addSquareRf input, switched to iLevel_vec.
% Updated 3/26/14 by DJ - fixed channel grid transpose issue, added outputs

% Set up
if ~exist('chansToPlot','var') || isempty(chansToPlot)
    chansToPlot = {'Fz';'Cz';'Pz';'Oz'};
end
if ~exist('iLevel_vec','var') || isempty(iLevel_vec)  
    iLevel_vec = repmat(length(R_cell{1}(1).responseFns),1,numel(R_cell));
end
if ~exist('addSquareRf','var') || isempty(addSquareRf)  
    addSquareRf = true;
end
% Get
nExp = numel(R_cell);
[RF,meanRF] = deal(cell(1,nExp));
for iExp=1:nExp
    iLevel = iLevel_vec(iExp);
    iEvents = iEvents_cell{iExp};
    R = R_cell{iExp};

    RF_cell = cell(1,numel(R));
    for i=1:numel(R)
        if iLevel>1 && addSquareRf
            RF_cell{i} = repmat(R(i).responseFns{1}(:,:,1),...
                [1,1,numel(iEvents)]) + R(i).responseFns{iLevel}(:,:,iEvents);
        else
            RF_cell{i} = R(i).responseFns{iLevel}(:,:,iEvents);
        end
    end

    RF{iExp} = cat(4,RF_cell{:});
    meanRF{iExp} = mean(RF{iExp},4);
    legendstr = R(1).regressor_events{iLevel}(iEvents);
    Cmap = GetSquaresEventColormap(legendstr);
    if iscell(R(1).tResponse)
        tResponse = R(1).tResponse{iLevel};
    else
        tResponse = R(1).tResponse;
    end
    chanlocs = R(1).EEG.chanlocs;
    
        figure(450+iExp);
        clf;

        PlotResponseFnsGrid(RF{iExp}, legendstr,tResponse,chanlocs,chansToPlot,Cmap);
        % PlotResponseFnsGrid(meanRF{iExp}(:,:,iEvents,:), legendstr(iEvents),tResponse,chanlocs,chanCol,Cmap(iEvents,:));
        nCols = size(chansToPlot,2);
        for i=1:numel(chansToPlot)
            subplot(size(chansToPlot,1),size(chansToPlot,2),i);
            iRow = floor((i-1)/nCols)+1;
            iCol = i-(iRow-1)*nCols;
            box on;
            ylabel(chansToPlot{iRow,iCol});
            title('');
            xlabel('');
            ylim([-4 4]);
        end
        xlabel('time (ms)');
        iDash = find(R(1).dataset=='-',1);
        prefix = R(1).dataset(1:iDash-1);
        MakeFigureTitle(sprintf('%s (%d subjects)',prefix,numel(R)));
end

% Resize and Cascade figures
ResizeFigures(450+(1:nExp),[],[1 nExp])