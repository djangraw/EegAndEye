% PlotContrastComparisons_script
%
% Created 9/17/14 by DJ.

%% Get info
experiment = 'sf3';
subject = 9;
old_suffix = 'GLMresults-Type-v3pt6-RampUp';
R = load(sprintf('%s-%d-%s',experiment,subject,old_suffix));
tResponse = R.tResponse{end};
chanlocs = R.EEG.chanlocs;
clear R;

%% PLOT

% experiments = 'sf';
experiments = {'sq','sf','sf3'};
% rules = {'D2vD0','D1vD0','D*vD0';'T0vD0','T1vD0','T2vD0','T+vD0'};
% rules = {'D2vD0','D1vD0','D*vD0','D+vD0'};%,'D-vD0'};
% rules = {'T0vD0','T1vD0','T2vD0','T+vD0'};%,'T-vD0'};
rules = 'T1vT0';

suffixes = 'Type-v3pt6-10fold';
% suffixes = {'Type-v3pt6-10fold','Type-v3pt6-Peak-10fold','Type-v3pt6-RampUp-10fold'};
multcorrect = 'fdr';
chansToPlot = {'FZ';'CZ';'PZ';'OZ'};
tBinCenters = 25:50:750;%[];
tBinWidth = 50;%[];

% Cmap = GetSquaresEventColormap({'pD_{2/3}','pD_{1/3}','pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pT_{+/3}'});
% Cmap = GetSquaresEventColormap({'pD_{2/3}','pD_{1/3}','pD^*_{-/3}','pD_{+/3}','pD_{-/3}'});
% Cmap = GetSquaresEventColormap({'pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pT_{+/3}','pT_{-/3}'});
Cmap = [];
smooth_std = 2.0;

[group_Z, group_RF, legendstr, titlestr] = CompileContrasts(experiments,suffixes,rules,multcorrect);
group_Z_smooth = SmoothData(group_Z,smooth_std,'full');
group_RF_smooth = SmoothData(group_RF,smooth_std,'full');

PlotResultsAsErps(group_RF_smooth,group_Z_smooth,chanlocs,tResponse,legendstr,Cmap,chansToPlot,titlestr);
PlotResultsAsScalpMaps(group_RF,group_Z,chanlocs,tResponse,tBinCenters,tBinWidth,legendstr,titlestr)


%% Make bar plot for SF3 experiment

experiments = 'sf3';
% rules = {'D2vD0','D1vD0','D*vD0','D+vD0','T0vD0','T1vD0','T2vD0','T+vD0'};
rules = {'D2vD0','D1vD0','T0vD0','T1vD0','T2vD0'};
suffixes = 'Type-v3pt6-10fold';
% Cmap = GetSquaresEventColormap({'pD_{2/3}','pD_{1/3}','pD^*_{-/3}','pD_{+/3}','pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pT_{+/3}'});
Cmap = GetSquaresEventColormap({'pD_{2/3}','pD_{1/3}','pT_{0/3}','pT_{1/3}','pT^*_{2/3}'});
multcorrect = 'fdr';

[group_Z, group_RF, legendstr, titlestr] = CompileContrasts(experiments,suffixes,rules,multcorrect);

% get component of max power in T2vD0 contrast
iComp = find(strcmp(rules,'T2vD0'));
[~,iMax] = max(sum(group_RF(:,:,iComp).^2,1));
compSM = group_RF(:,iMax,iComp);
% normalize
compSM = compSM/sqrt(sum(compSM.^2));

compTC = zeros(1,size(group_RF,2),size(group_RF,3));
compZ = zeros(1,size(group_RF,2),size(group_RF,3));
for i=1:size(group_RF,3)
    compTC(1,:,i) = compSM'*group_RF(:,:,i);     
    compZ(1,:,i) = compSM'*group_Z(:,:,i)/sqrt(sum(compSM.^2));
end

compTC_smooth = SmoothData(compTC,smooth_std,'full');
compZ_smooth = SmoothData(compZ,smooth_std,'full');
PlotResultsAsErps(compTC_smooth,compZ_smooth,chanlocs(1),tResponse,legendstr,Cmap,{chanlocs(1).labels},titlestr);

set(gcf,'Position',[1 1133 1147 373]);

%% Add window markers
iWin = (iMax-2):(iMax+2); % 50ms window
PlotVerticalLines(tResponse(iWin([1 end])),'k--');


%% 
figure(214); clf;
axis; hold on;


for i=1:size(compTC,3)
    bar(i,squeeze(mean(compTC(1,iWin,i),2)),'facecolor',Cmap(i,:));
end
set(gca,'xtick',1:length(rules),'xticklabel',rules);
set(gcf,'Position',[1 800 1147 373]);



    





