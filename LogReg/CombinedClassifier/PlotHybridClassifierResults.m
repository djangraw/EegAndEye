function PlotHybridClassifierResults(classifierNames,chanlocs,varargin)

% PlotHybridClassifierResults(classifierNames,chanlocs,results1,results2,...)
%
% INPUTS:
% -classifierNames is a cell array of strings.
% -chanlocs is an n-element cell array of EEG.chanlocs structs, where n is
% the number of subjects.
% -results1, etc. are n-element arrays of results structs from 
% ClassifyWithEegAndDwellTime.m. Each should contain information about the 
% classification, with fields w,v,fwdModel,Az.
%
% OUTPUTS:
% -plots a figure for each classifier with the fwd models, spatial weights, etc.
%
% Created 7/11/13 by DJ.

%% Declare constants
eeg_bintimes = 100:100:900; % in ms
eeg_binwidth = 100; % in ms
eeg_binctr = eeg_bintimes + eeg_binwidth/2;
nClassifiers = length(varargin);
nSubjects = length(varargin{1});
ps_bintimes = 0:500:2500;
ps_binwidth = 500;
ps_binctr = ps_bintimes + ps_binwidth/2;
nPSbins = numel(ps_bintimes);
% Old way - put DT to right of PS
dt_binctr = ps_binctr(end)+ps_binctr(end)-ps_binctr(end-1);
stderr_DT = 0;
% New way - use DT averages as x position
load('TEMP_DwellTimes.mat');   
dwelltimes = nan(1,nSubjects);
for j=1:nSubjects
    dwelltimes(j) = mean(DT{j});
end        
dt_binctr = mean(dwelltimes);
stderr_DT = std(dwelltimes)/sqrt(nSubjects);




%% Get Channel Locations
if isnumeric(chanlocs)
    subjects = chanlocs;
    disp('Getting channel locations')
    chanlocs = cell(1,nSubjects);
    for i=1:nSubjects
        foo = load(sprintf('ALLEEG%d_NoeogEpochedIcaCropped.mat',subjects(i)));
        chanlocs{i} = foo.ALLEEG(1).chanlocs;
    end
end

% Reconcile different channels being used (use subset common to all subjs)
channies = [chanlocs{:}];
chans = unique({channies.labels});
for i=1:nSubjects
    chans = chans(ismember(chans,{chanlocs{i}.labels}));
end
chanlocs_all = chanlocs{1}(ismember({chanlocs{1}.labels},chans));

%% Set plot parameters
fontname = 'Futura'; % Futura, Helvetica, Gotham are all good
fontsize = 15;
linewidth = 2;
markersize = 20;

%% Main loop
Az_all = nan(nClassifiers,nSubjects);
for iClassifier = 1:nClassifiers
    R = varargin{iClassifier};
    % extract
    w_raw = cat(5,R.w);    
    v_raw = cat(4,R.v);
    Az_raw = cat(2,R.Az);
    % mean & stderr across folds & subjects
    w = mean(mean(mean(w_raw,3),4),5);
    v = mean(mean(v_raw,3),4);
    stderr_w = std(mean(mean(w_raw,3),4),[],5)/sqrt(nSubjects);
    stderr_v = std(mean(v_raw,3),[],4)/sqrt(nSubjects);
    
    % fwd models
    fwds = zeros(numel(chans), size(R(1).fwdModel,2),nSubjects);    
    for i=1:nSubjects
       fwds(:,:,i) = mean(R(i).fwdModel(ismember({chanlocs{i}.labels},chans),:,:),3);
    end
    
    % separate PS & DT
    nEegBins = size(w,2);
    v_EEG = v(1:nEegBins);
    stderr_v_EEG = stderr_v(1:nEegBins);
    v(1:nEegBins) = [];
    stderr_v(1:nEegBins) = [];
    nOtherBins = size(v,2);
    if nOtherBins==0
        [v_DT, v_PS, stderr_v_DT, stderr_v_PS] = deal([]);        
    elseif nOtherBins==1;
        v_DT = v;        
        v_PS = [];
        stderr_v_DT = stderr_v;
        stderr_v_PS = [];
    elseif nOtherBins==nPSbins
        v_PS = v;
        v_DT = [];        
        stderr_v_PS = stderr_v;
        stderr_v_DT = [];
    elseif nOtherBins==nPSbins+1
        v_PS = v(1:end-1);
        v_DT = v(end);
        stderr_v_PS = stderr_v(1:end-1);
        stderr_v_DT = stderr_v(end);
    else
        error('Check hard-coded nPSbins!')
    end    
    
    % Set up plot
    figure(300+iClassifier); clf;
    
    %% Plot scalp maps
    % cmax = max(max(abs(mean(fwds,3))));            
    nBins = numel(eeg_bintimes);
    nCols = ceil(nBins/2);
    cmax = 10;
    for i=1:size(fwds,2)
        subplot(2,nCols,rem(i-1,2)*nCols+floor((i-1)/2)+1,'fontname',fontname,'fontsize',fontsize)
        topoplot(mean(fwds(:,i,:),3),chanlocs_all,'maplimits',[-cmax cmax],'plotrad',0.5,'electrodes','off'); % don't include mullet
    %     topoplot(mean(fwds(:,i,:),3),chanlocs_all,'maplimits',[-cmax cmax],'conv','on'); % do include mullet
%         title(sprintf('%d-%dms',eeg_bintimes(i),eeg_bintimes(i)+eeg_binwidth));
        title(sprintf('%d ms',eeg_binctr(i)));
    end
    % Annotate plot
    axes('Position',[.85 .55 .1 .24],'CLim',[-cmax cmax],'visible','off','fontname',fontname,'fontsize',fontsize);
    colorbar('fontname',fontname,'fontsize',fontsize)
    % Annotate figure
    h = MakeFigureTitle('Forward Models (Mean Across Subjects) (uV)');
    set(h,'fontname',fontname,'fontsize',fontsize,'fontweight','normal')

    %% Plot temporal weights
    figure(400+iClassifier); clf;
    set(gca,'fontname',fontname,'fontsize',fontsize); hold on;
%     subplot(2,1,2,'fontname',fontname,'fontsize',fontsize); hold on;

    % Plot temporal weights
    legendstr = {};
    if ~isempty(v_EEG)
        errorbar(eeg_binctr(1:length(v_EEG)),v_EEG,stderr_v_EEG,...
            'b.-','linewidth',linewidth,'markersize',markersize);
        legendstr = [legendstr {'EEG'}];
    end
    if ~isempty(v_PS)
        errorbar(ps_binctr(1:length(v_PS)),v_PS,stderr_v_PS,...
            '.-','linewidth',linewidth,'markersize',markersize,...
            'color',[0 0.5 0]);
        legendstr = [legendstr {'Pupil Dilation'}];
    end
    if ~isempty(v_DT)
        errorbar(dt_binctr(1:length(v_DT)),v_DT,stderr_v_DT,...
            'r.-','linewidth',linewidth,'markersize',markersize);    
        hh = herrorbar(dt_binctr(1:length(v_DT)),v_DT,stderr_DT,'r.-');
        set(hh,'linewidth',linewidth,'markersize',markersize);
        legendstr = [legendstr {'Dwell Time'}];
    end
    
    % Set xlims
    if ~isempty([v_DT v_PS])
        xmax = max([dt_binctr ps_binctr])+ps_binwidth/2;
    elseif ~isempty(v_PS)
        xmax = max(ps_binctr)+ps_binwidth/2;
    else
        xmax = max(eeg_binctr)+eeg_binwidth/2;
    end
        
    set(gca,'xlim',[0 xmax]);

    % Annotate plot
    set(gca,'xgrid','on','box','on')
    ylabel('Second-level Feature Weights')
    xlabel('time from saccade offset to bin center (ms)')    
    legend(legendstr,'Location','NorthEast')
    hold on
    plot(get(gca,'xlim'),[0 0],'k--','LineWidth',linewidth);
    
    % Record
    Az_all(iClassifier,:) = Az_raw;

end

%% Plot Az's
figure(300+nClassifiers+1);
clf; hold on;
set(gca,'xtick',1:nSubjects,'ytick',0.3:.1:1,'box','on','fontname',fontname,'fontsize',fontsize)
% Plot
PlotUniqueLines(1:nSubjects, Az_all', '.', linewidth, markersize)
plot([0 nSubjects+1],[0.5 0.5],'k--','LineWidth',linewidth)
% Annotate
xlabel('subject')
ylabel('Area Under ROC Curve')
title('Classifier Performance')
ylim([0.3 1])
legend(classifierNames{:})