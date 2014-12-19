function [density,density_all,plotX,plotY] = NEDE_MakeEyeHeatmap(y,nBins,clim)

% Plot heatmap of eye position for each session and across sessions.
%
% [density,density_all,plotX,plotY] = NEDE_MakeEyeHeatmap(y,nBins,clim)
%
% INPUTS:
% - y is an n-element vector of NEDE structs.
% - nBins is a scalar indicating the number of bins into which you want to
% separate the density plots. More bins = finer resolution. [default = 128]
% - clim is a 2-element vector indicating the min and max you would like to
% make the color limits in each plot (in units of the proportion of all 
% time points). [default = [0 1e-4]]
%
% OUTPUTS:
% - density is an n-element vector of cells, in which density{i} is an
% nBins x nBins matrix of the eye position density in session i.
% - density_all is an nBins x nBins matrix of the eye position density
% across all sessions.
% - plotX and plotY are 1 x nBins vectors of the bin centers being plotted.
%
% Created 11/7/14 by DJ.

% Declare defaults
if ~exist('nBins','var') || isempty(nBins)
    nBins = 128;
end
if ~exist('clim','var') || isempty(clim)
    clim = [0 1e-4];
end

% --- Set up
% Declare constants
nPlots = numel(y)+1;
nRows = ceil(sqrt(nPlots));
nCols = ceil(nPlots/nRows);
MIN_XY = [0 0];
MAX_XY = [y(1).params.screen.width, y(1).params.screen.height];
% Set up loop
data = cell(1,numel(y));
density = cell(1,numel(y));
[meanpos, medianpos] = deal(nan(numel(y),2));

% --- Main loop
for i=1:numel(y)
    % Calculate density
    data{i} = y(i).events.fixupdate.position;
    [~,density{i},X,Y] = kde2d(data{i},nBins,MIN_XY,MAX_XY);
    subplot(nRows,nCols,i);
    imagesc(X(1,:),Y(:,1),density{i});
    % Get mean and median eye pos
    meanpos(i,:) = nanmean(data{i},1);
    medianpos(i,:) = nanmedian(data{i},1);    
    % superimpose mean eye pos
    hold on;
%     plot(meanpos(i,1),meanpos(i,2),'w+');
    plot(medianpos(i,1),medianpos(i,2),'r+','markersize',10);
    axis([0 MAX_XY(1) 0 MAX_XY(2)]); % put limits back
%     axis equal
    % annotate plot
    colorbar
    set(gca,'clim',clim,'ydir','reverse');
    title(show_symbols(x.params.EDF_filename(1:end-4)));    
end

%% Average
subplot(nRows,nCols,numel(sessions)+1); cla;
[~,density_all,X,Y] = kde2d(cat(1,data{:}),nBins,MIN_XY,MAX_XY);

% Plot averagre
imagesc(X(1,:),Y(:,1),density_all);
% annotate plot
colorbar
set(gca,'clim',clim);
prefix = x.params.EDF_filename(1:find(x.params.EDF_filename=='-',1)-1);
title(show_symbols(sprintf('%s-%d, all %d sessions',prefix,x.params.subject,numel(y))));
% create output
plotX = X(1,:);
plotY = Y(:,1)';