function PlotResultsAsScalpMaps(group_RF,group_Z,chanlocs,tResponse,tBinCenters,tBinWidth,legendstr,titlestr)

% PlotResultsAsScalpMaps(group_RF,group_Z,chanlocs,tResponse,tBinCenters,tBinWidth,legendstr,titlestr)
%
% Created 9/16/14 by DJ based on GetContrastResults_fast.

if ~exist('tBinCenters','var') || isempty(tBinCenters)
    tBinCenters = 37.5:75:1000;%750;
end
if ~exist('tBinWidth','var') || isempty(tBinWidth)
    tBinWidth = 75;
end

% Plot RFs
if ~isempty(group_RF)
    figure(222); clf;
    scalpmaps = GetScalpMaps(group_RF,tResponse,tBinCenters,tBinWidth);
    PlotScalpMaps(scalpmaps,chanlocs,[],tBinCenters-tBinWidth/2,legendstr);
    % Annotate figure
    MakeFigureTitle(sprintf('%s, Group RFs',titlestr));
    nRows = size(group_RF,3);
    set(gcf,'Position',[825 1506-nRows*80 1261 nRows*80]);
end

% Plot Z scores
cthresh = 1.96; % z score for 2-tailed p=0.05

if ~isempty(group_Z)
    figure(223); clf;
    scalpmaps = GetScalpMaps(group_Z,tResponse,tBinCenters,tBinWidth,cthresh);
    PlotScalpMaps(scalpmaps,chanlocs,[],tBinCenters-tBinWidth/2,legendstr);
    % Annotate figure
    MakeFigureTitle(sprintf('%s, Group Z scores',titlestr));
    set(gcf,'Position',[825 1506-2*nRows*80 1261 nRows*80]);
end
