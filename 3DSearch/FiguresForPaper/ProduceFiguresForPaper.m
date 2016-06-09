% ProduceFiguresForPaper.m
%
% Created 6/7/13 by DJ.
% Updated several times...
% Updated 12/31/13 by DJ - switched to TEMP_HybridClassifier_Results,
%   added response to reviewers figures.

%% SET UP
% Set plot parameters
fontname = 'Futura'; % Futura, Helvetica, Gotham, Georgia are all good
fontsize = 12;
linewidth = 2;
markersize = 20;

%% Declare constants
subjects = [22:30 32];
folders = {'2013-03-27-Pilot-S22', '2013-03-28-Pilot-S23', ...
    '2013-04-29-Pilot-S24', '2013-05-01-Pilot-S25', '2013-05-02-Pilot-S26',...
    '2013-05-03-Pilot-S27', '2013-05-07-Pilot-S28', '2013-05-10-Pilot-S29',...
    '2013-05-15-Pilot-S30', '2013-06-06-Pilot-S32'};
sessions_cell = {2:14, [3 6:17], 1:15, 1:15, 1:15, 1:15, 1:15, 1:15, [1:10 12:15], 2:16};
offsets = [-12 -4 48 60 68 68 92 112 -32 88];
fracSensitivity = 0.9999999; % used to calculate nSensitivity
homedir = '/Users/dave/Documents/Data/3DSearch';
homedir = '/Users/jangrawdc/Documents/LiincData/3DSearch';
cvmode = '10fold';
usegridconstraints = true;
levelname = 'GridHuge.jpg';

%% Load Data
% load TEMP_EegAndDwell_Results.mat
load TEMP_HybridClassifier_Results.mat

dwell = cell(1,numel(subjects));
isToTarget = cell(1,numel(subjects));
for i=1:numel(subjects)
    fprintf('--- Subject %d/%d: ---\n',i,numel(subjects));
    [dwell{i} isToTarget{i}] = GetDwellTimes(subjects(i),sessions_cell{i});
    dwellz = [dwell{i}{:}];
    isTargetz = [isToTarget{i}{:}];
    y_dwell{i} = dwellz(~isnan(dwellz));
    truth_dwell{i} = isTargetz(~isnan(dwellz));
end
all_dwell = [y_dwell{:}];
all_truth = [truth_dwell{:}];

%% Dwell Time figure

nSubjects = numel(subjects);
targs_at_dwell = zeros(numel(subjects),max(all_dwell));
dists_at_dwell = zeros(numel(subjects),max(all_dwell));
for iSubj = 1:numel(subjects)
    for i=1:max(all_dwell)
        targs_at_dwell(iSubj,i) = mean(y_dwell{iSubj}(truth_dwell{iSubj})>i);
        dists_at_dwell(iSubj,i) = mean(y_dwell{iSubj}(~truth_dwell{iSubj})>i);
    end    
end

figure(247); cla; hold on;
% plot([targs_at_dwell; dists_at_dwell]','LineWidth',2);
% plot([targs_at_dwell-dists_at_dwell]');
xDwell = 1:max(all_dwell);
plot(xDwell,mean(dists_at_dwell),'b','LineWidth',2);
plot(xDwell,mean(targs_at_dwell),'r','LineWidth',2);
ErrorPatch(xDwell,mean(dists_at_dwell),std(dists_at_dwell)/sqrt(nSubjects),'b','b');
ErrorPatch(xDwell,mean(targs_at_dwell),std(targs_at_dwell)/sqrt(nSubjects),'r','r');

% Annotate plot
set(gca,'box','on','fontname',fontname,'fontsize',fontsize);
xlabel('t (ms)');
ylabel(sprintf('Fraction of Trials with Dwell Times > t\n(Mean +/- std err across subjects)'))
legend('Distractors','Targets');

%% Pupil Size Figure
figure;
set(gcf,'Position',[1000 992 733 486]);
% subplot(2,1,1); 
clf; hold on;
set(gca,'box','on','fontname',fontname,'fontsize',fontsize)
mtarg = mean(ps_median_targ_all,1);
mdist = mean(ps_median_dist_all,1);
mdiff = mean(ps_median_targ_all - ps_median_dist_all,1);
starg = std(ps_median_targ_all,[],1)/sqrt(numel(subjects));
sdist = std(ps_median_dist_all,[],1)/sqrt(numel(subjects));
sdiff = std(ps_median_targ_all - ps_median_dist_all,[],1)/sqrt(numel(subjects));


plot(t_epoch,mdist,'b','linewidth',2);
plot(t_epoch,mtarg,'r','linewidth',2);
% plot(t_epoch,mdiff,'g','linewidth',2);
ErrorPatch(t_epoch,mdist,sdist,'b','b');
ErrorPatch(t_epoch,mtarg,starg,'r','r');
% ErrorPatch(t_epoch,mdiff,sdiff,'g','g');
plot([t_epoch(1) t_epoch(end)],[0 0],'k-');
plot([0 0],get(gca,'ylim'),'k-');
set(colorbar,'Visible','off') % dummy colorbar to make axes line up
% legend('Target','Distractor','Difference','Location','NorthWest');
legend('Distractors','Targets','Location','NorthEast');
ylabel(sprintf('Median pupil area change\n(%% of subj mean)\nmean +/- stderr across subjects'))
xlabel('Time from first saccade to object (ms)');

%% Plot ERPs
gridOfChans = {'F5' 'FZ' 'F6'; 'C5' 'CZ' 'C6'; 'P5' 'PZ' 'P6'};
[all_dist,all_targ,times,chanlocs_all] = TEMP_Make3dsErps(subjects,sessions_cell,offsets,gridOfChans);

%%

gridOfChans = {'FZ';'CZ';'PZ';'O1'};
% Re-subtract baseline
% tBase = [0 100];
% iBaseTimes = find(times>tBase(1) & times<tBase(2));
% for i=1:size(all_dist,4)    
%     debased_dist(:,:,:,i) = all_dist(:,:,:,i) - repmat(mean(all_dist(:,iBaseTimes,:,i),2),[1 length(times)]);
%     debased_targ(:,:,:,i) = all_targ(:,:,:,i) - repmat(mean(all_targ(:,iBaseTimes,:,i),2),[1 length(times)]);
% end

% Re-plot
% PlotResponseFnsGrid(cat(3,debased_dist,debased_targ),{'Distractors', 'Targets'},times,chanlocs_all,gridOfChans);
% PlotResponseFnsGrid(cat(3,all_dist,all_targ),{'Distractors', 'Targets'},times,chanlocs_all,gridOfChans);
RFdiff = SmoothData(all_targ-all_dist,6,'full');
PlotResponseFnsGrid(RFdiff,{'Targets - Distractors'},times,chanlocs_all,gridOfChans);
% Annotate
set(gcf,'Position',[0 623 704 882]);
for i=1:4
    subplot(4,1,i);
    title('');    
    ylabel(gridOfChans{i});
    set(gca,'xtick',-200:100:1000);
    if i<3
        set(gca,'xticklabel',{});
        xlabel('');
    end
    xlim([-200 1000]);
%     ylim([-2 5])
    ylim([-1 2])
end

%% Plot scalp maps
mean_dist = mean(all_dist,4);
mean_targ = mean(all_targ,4);
mean_diff = mean(all_targ-all_dist,4);
figure;
tPlots = -50:100:1000;
binWidth = 100;
responses = GetScalpMaps(cat(3,mean_dist,mean_targ,mean_diff),times,tPlots,binWidth);
PlotScalpMaps(responses,chanlocs_all,[],tPlots,{'Dist','Targ','Diff'});

%% Get chanlocs
chanlocs = cell(1,numel(subjects));
for i=1:numel(subjects)
    fprintf('--- Subject %d/%d ---\n',i,numel(subjects));
    load(sprintf('ALLEEG%d_NoeogEpochedIcaCropped.mat',subjects(i))); % ALLEEG
    chanlocs{i} = ALLEEG(1).chanlocs;
end
disp('Done!')

%% Classifier fwd models/weights figure

% Declare which data to use
R = R_eegpsdt;
classifierNames = {'Hybrid'};

PlotHybridClassifierResults(classifierNames,chanlocs,R)
figure(301);
set(gcf,'Position',[0 810 1220 800]);
subplot(2,1,2);
set(gca,'ylim',[-0.2 0.9]);
set(gca,'ytick',[-0.2:0.2:0.8]);



%% PLOT AZS

Az_all = [Az_eeg; Az_ps; Az_dt; Az_eegpsdt];
classifierNames = {'EEG','Pupil Dilation','Dwell Time','Hybrid','Chance'};

% Az_all = [Az_ps; [R_ps_late.Az]; Az_eegpsdt; [R_eegpsdt_late.Az]];
% classifierNames = {'PD','PD (1-3s only)','Hybrid','Hybrid (1-3s only)','Chance'};
% Az_all = [Az_eeg; [R_ps_late.Az]; Az_dt; [R_eegpsdt_late.Az]];
% classifierNames = {'EEG','Pupil Dilation (1-3s only)','Dwell Time','Hybrid (1-3s only)','Chance'};

% Az_all = [Az_eegps; Az_eegdt; Az_psdt; Az_eegpsdt];
% classifierNames = {'EEG+PD','EEG+DT','PD+DT','Hybrid','Chance'};
% Az_all = [Az_eeg; Az_ps; Az_dt; Az_eegpsdt; Az_eegpsdt_cv];
% classifierNames = {'EEG','Pupil Area','Dwell Time','Hybrid','Hybrid + CV','Chance'};
nSubjects = numel(subjects);
[~,order] = sort(Az_eeg,'descend');
Az_all = Az_all(:,order);

figure(250); clf; hold on;
% set(gcf,'Position',[2   912   836   572]);
set(gcf,'Position',[2        1179         443         305]);
set(gca,'xtick',1:nSubjects,'ytick',0.3:.1:1,'box','on','fontname',fontname,'fontsize',fontsize)
% Plot
% PlotUniqueLines(1:nSubjects, Az_all', '.', linewidth, markersize)
plot(1:nSubjects, Az_all', '.-','linewidth',linewidth,'markersize',markersize)
plot([0 nSubjects+1],[0.5 0.5],'k--','LineWidth',linewidth)
% Annotate
xlabel('Subject')
ylabel('AUC')
% xlabel('Subject (sorted by EEG Az score)')
% ylabel('Area Under ROC Curve')
% title('Classifier Performance')
ylim([0.3 1])
legend(classifierNames{:})
xlim([0 nSubjects+1])


%% TSP Route Figure (Get)
% Declare which data you want to use
% R = R_both_noeog_eica20_nested;
R = R_eegpsdt;
TagTours = TagTours_eegpsdt;
plotInfos = plotInfos_eegpsdt;
iSubj = order(7);

% % Get Predictions
% load(sprintf('ALLEEG%d_NoeogEpochedIcaCropped.mat',subjects(iSubj))); % ALLEEG
% y = loadBehaviorData(subjects(iSubj),sessions_cell{iSubj},'3DS');
% iObjects_eeg_pt{iSubj} = GetEegPredictedTargets(subjects(iSubj),sessions_cell{iSubj},ALLEEG,R(iSubj).y,R(iSubj).truth);

% Get Stats and TSP route
if exist('TagTours','var') && exist('plotInfos','var')
    TagTour = TagTours{iSubj};
    plotInfo = plotInfos{iSubj};
else
    [~,TagTour,plotInfo] = CalculateImprovement(subjects(iSubj),sessions_cell{iSubj},levelname,iObjects_eeg_pt{iSubj},usegridconstraints,fracSensitivity);
end
%% (Plot)
TagTour_grid = [];
for i=1:size(TagTour,1)-1
    TagTour_grid = [TagTour_grid; TagTour(i,:)];
    if TagTour(i,2)<TagTour(i+1,2)
        TagTour_grid = [TagTour_grid; TagTour(i,1), TagTour(i,2)+10];
        TagTour_grid = [TagTour_grid; TagTour(i+1,1), TagTour(i,2)+10];
    elseif TagTour(i,2)>TagTour(i+1,2)
        TagTour_grid = [TagTour_grid; TagTour(i,1), TagTour(i,2)-10];
        TagTour_grid = [TagTour_grid; TagTour(i+1,1), TagTour(i,2)-10];
    elseif size(TagTour_grid,1)==1 || TagTour_grid(end,2)>TagTour_grid(end-1,2)
        TagTour_grid = [TagTour_grid; TagTour(i,1), TagTour(i,2)+10];
        TagTour_grid = [TagTour_grid; TagTour(i+1,1), TagTour(i,2)+10];
    else
        TagTour_grid = [TagTour_grid; TagTour(i,1), TagTour(i,2)-10];
        TagTour_grid = [TagTour_grid; TagTour(i+1,1), TagTour(i,2)-10];
    end
%     TagTour_grid = [TagTour_grid; TagTour(i,1), TagTour(i+1,2)];
end
TagTour_grid = [TagTour_grid; TagTour(end,:)];
rColor = repmat(0.75,1,3);
figure(248); clf;
% set(gcf,'Position',[1 821 2560 663]);
set(gcf,'Position',[1         853        1076         631]);

for i=1:6
    subplot(3,2,i);
    axis equal
%     % view small window
%     xlim([500 1300]);
%     ylim([-20 180]);
    % view ~1/4 of enviro
    % For subj 22
    xlim([-9 427]);
    ylim([5 155]);
    % For subj 25
%     xlim([590 1030]);
%     xlim([562 998]);
%     ylim([5 155]);
    set(gca,'xtick',[],'ytick',[],'box','on','fontname',fontname,'fontsize',fontsize);
end

% Get info
wasseen = ismember(plotInfo.objlocs,cat(1,plotInfo.cameraPath{:}),'Rows');
dotsize = (length(plotInfo.tagRank)+1-plotInfo.tagRank)*500/length(plotInfo.tagScore);
% dotsize = plotInfo.tagScore*500/length(plotInfo.tagScore);
isintsp = ismember(plotInfo.objlocs,TagTour,'Rows');

subplot(3,2,1); cla; hold on;
% title('A. hBCI Classification')
% % Plot line for camera path
% for i=1:numel(plotInfo.cameraPath)
%     plot(plotInfo.cameraPath{i}(:,1),plotInfo.cameraPath{i}(:,2),'c--');
% end
% Plot patches for camera path
for i=1:numel(plotInfo.cameraPath)
    thispath = plotInfo.cameraPath{i};
    cols = unique(thispath(:,1));
    for j=1:numel(cols)
        rows = thispath(thispath(:,1)==cols(j),2);
        rectangle('Position',[cols(j)-7.5, min(rows)-10, 15, range(rows)+20],'FaceColor',rColor,'EdgeColor',rColor);
    end
end

scatter(plotInfo.objlocs(:,1),plotInfo.objlocs(:,2),'b.');
scatter(plotInfo.objlocs(plotInfo.iTargets,1),plotInfo.objlocs(plotInfo.iTargets,2),'r.');
% scatter(plotInfo.objlocs(~wasseen,1),plotInfo.objlocs(~wasseen,2),'k.')
plot(plotInfo.objlocs(plotInfo.iEegTargets,1),plotInfo.objlocs(plotInfo.iEegTargets,2),'mo','linewidth',linewidth,'markersize',8);

subplot(3,2,3); cla; hold on;
% title('B. CV Extrapolation')
% Plot patches for camera path
for i=1:numel(plotInfo.cameraPath)
    thispath = plotInfo.cameraPath{i};
    cols = unique(thispath(:,1));
    for j=1:numel(cols)
        rows = thispath(thispath(:,1)==cols(j),2);
        rectangle('Position',[cols(j)-7.5, min(rows)-10, 15, range(rows)+20],'FaceColor',rColor,'EdgeColor',rColor);
    end
end
scatter(plotInfo.objlocs(:,1),plotInfo.objlocs(:,2),dotsize,'b.');
scatter(plotInfo.objlocs(plotInfo.iTargets,1),plotInfo.objlocs(plotInfo.iTargets,2),dotsize(plotInfo.iTargets),'r.');
% scatter(plotInfo.objlocs(~wasseen,1),plotInfo.objlocs(~wasseen,2), dotsize(~wasseen),'k.')
plot(plotInfo.objlocs(plotInfo.iTagTargets,1),plotInfo.objlocs(plotInfo.iTagTargets,2),'s','color',[0 .5 0],'linewidth',linewidth,'markersize',10);

subplot(3,2,5); cla; hold on;
% title('C. Route Planning')
% Plot patches for camera path
for i=1:numel(plotInfo.cameraPath)
    thispath = plotInfo.cameraPath{i};
    cols = unique(thispath(:,1));
    for j=1:numel(cols)
        rows = thispath(thispath(:,1)==cols(j),2);
        rectangle('Position',[cols(j)-7.5, min(rows)-10, 15, range(rows)+20],'FaceColor',rColor,'EdgeColor',rColor);
    end
end
% dotsize = plotInfo.tagScore*500/length(plotInfo.tagScore);
scatter(plotInfo.objlocs(:,1),plotInfo.objlocs(:,2),dotsize,'b.');
scatter(plotInfo.objlocs(plotInfo.iTargets,1),plotInfo.objlocs(plotInfo.iTargets,2),dotsize(plotInfo.iTargets),'r.');
% scatter(plotInfo.objlocs(~wasseen & ~isintsp,1),plotInfo.objlocs(~wasseen & ~isintsp,2),dotsize(~wasseen & ~isintsp),'k.')
plot(plotInfo.objlocs(plotInfo.iTagTargets,1),plotInfo.objlocs(plotInfo.iTagTargets,2),'s','color',[0 .5 0],'linewidth',linewidth,'markersize',10);
% plot(TagTour(:,1),TagTour(:,2),'k--');
plot(TagTour_grid(:,1),TagTour_grid(:,2),'k--','linewidth',linewidth);

subplot(3,2,4); cla; hold on;
title('Truth Data')
% Plot patches for camera path
for i=1:numel(plotInfo.cameraPath)
    thispath = plotInfo.cameraPath{i};
    cols = unique(thispath(:,1));
    for j=1:numel(cols)
        rows = thispath(thispath(:,1)==cols(j),2);
        rectangle('Position',[cols(j)-7.5, min(rows)-10, 15, range(rows)+20],'FaceColor',rColor,'EdgeColor',rColor);
    end
end
% dotsize = plotInfo.tagScore*500/length(plotInfo.tagScore);
% scatter(plotInfo.objlocs(:,1),plotInfo.objlocs(:,2),dotsize,'k.');
scatter(plotInfo.objlocs(:,1),plotInfo.objlocs(:,2),dotsize,'b.');
scatter(plotInfo.objlocs(plotInfo.iTargets,1),plotInfo.objlocs(plotInfo.iTargets,2),dotsize(plotInfo.iTargets),'r.');
plot(plotInfo.objlocs(plotInfo.iEegTargets,1),plotInfo.objlocs(plotInfo.iEegTargets,2),'mo','linewidth',linewidth,'markersize',8);
plot(plotInfo.objlocs(plotInfo.iTagTargets,1),plotInfo.objlocs(plotInfo.iTagTargets,2),'s','color',[0 .5 0],'linewidth',linewidth,'markersize',10);
% plot(plotInfo.cameraPath{1}(:,1),plotInfo.cameraPath{1}(:,2),'c--');
patch([-1 -1 0],[-1 0 0],rColor,'EdgeColor',rColor); % dummy patch for legend
% plot(TagTour(:,1),TagTour(:,2),'k--');
plot(TagTour_grid(:,1),TagTour_grid(:,2),'k--','linewidth',linewidth);
% for i=1:numel(plotInfo.cameraPath)
%     plot(plotInfo.cameraPath{i}(:,1),plotInfo.cameraPath{i}(:,2),'c--');
% end
% Make legend
legend('Distractor','Target','hBCI predicted tartget','CV predicted target','Explored Areas','Traveling Salesman Route');
% legend('Unseen object','Distractor','Target','hBCI predicted tartget','CV predicted target','Explored Areas','Traveling Salesman Route');

%% SET CHANCE LEVELS
chance.precision = [25 25];
chance.pctFound = [20/56*.25*100 25];
chance.pctDistance = 100;

%% Get best performance metrics
for i=1:nSubjects
    [stats_best(i),TagTours_best{i},plotInfos_best{i}] = FindBestRoute(subjects(i),sessions_cell{i},levelname,usegridconstraints);
end

best.precision = [100 100];
best.pctFound = [20/56*100, 100];
best.pctDistance = mean([stats_best.pctDistance]);

%% System Performance Figure
% NOTE: LEGEND MUST BE MOVED BY HAND TO THE LEFT OF FIRST PLOT!

% Declare which stats to use
stats = stats_eegpsdt;
% stats = stats_dt;
[~,order] = sort(Az_eeg,'descend');

subjectstrings = cell(1,numel(subjects));
for i=1:numel(subjects)
    subjectstrings{i} = sprintf('S%d',i);
end

% Make plots
% figure(251); clf;
figure(254); clf;
% set(gcf,'Position',[1 540 1730 944]);
set(gcf,'Position',[3   912   747   572]);
subplot(2,2,1); cla; hold on;
PlotUniqueLines([1 2], cat(1,stats(order).precision)','.',linewidth,markersize,nSubjects);
plot(chance.precision,'k.:','linewidth',linewidth,'markersize',markersize);
plot(best.precision,'k.--','linewidth',linewidth,'markersize',markersize);
set(gca,'xtick',[1 2],'xticklabel',{'hBCI','hBCI+CV'},'box','on','fontname',fontname,'fontsize',fontsize)
xlim([0.5 2.5]);
ylim([0 100]);
ylabel('% Precision');
xlabel('Target Prediction System');
% title('Precision of Predicted Targets');
legend([subjectstrings, {'Chance', 'Perfect'}],'Location','SouthEast')

% figure(252); clf;
subplot(2,2,2); cla; hold on;
PlotUniqueLines([1,2],cat(1,stats(order).pctFound)','.',linewidth,markersize,nSubjects);
plot(chance.pctFound,'k.:','linewidth',linewidth,'markersize',markersize);
plot(best.pctFound,'k.--','linewidth',linewidth,'markersize',markersize);
% set(gca,'xtick',[1 2],'xticklabel',{'Training','TAG'},'box','on','fontname',fontname,'fontsize',fontsize)
set(gca,'xtick',[1 2],'xticklabel',{'hBCI','hBCI+CV'},'box','on','fontname',fontname,'fontsize',fontsize)
xlim([0.5 2.5]);
ylim([0 100]);
xlabel('Target Prediction System');
ylabel('% Targets Identified')
% title('Percentage of True Targets Identified');
% legend([subjectstrings, {'Chance'}],'Location','SouthEast')

% figure(253); clf;
subplot(2,2,3); cla; hold on;
pctFound = cat(1,stats(order).pctFound);
PlotUniqueLines([nan(numel(stats),1), cat(1,stats(order).pctDistance)]', [zeros(numel(stats),1), pctFound(:,2)]','.',linewidth,markersize,nSubjects);
plot([0 chance.pctDistance],[0 100], 'k:','linewidth',linewidth,'markersize',markersize);
plot([0 best.pctDistance],[0 100], 'k--','linewidth',linewidth,'markersize',markersize);
set(gca,'box','on','fontname',fontname,'fontsize',fontsize)
xlim([0 100]);
ylim([0 100]);
ylabel('% Targets Viewed')
xlabel('% Distance Traveled')
% title('Search Efficiency with Full System');
% legend([subjectstrings, {'View All'}],'Location','SouthEast')


%% SUPPLEMENTARY System Performance Figure
% NOTE: LEGEND MUST BE MOVED BY HAND.

% Declare which stats to use
% stats = stats_eegpsdt;
stats_cell = {stats_eeg, stats_ps, stats_dt, stats_eegpsdt};
stats_names = {'EEG','Pupil Dilation','Dwell Time','Hybrid'};
[~,order] = sort(Az_eeg,'descend');

subjectstrings = cell(1,numel(subjects));
for i=1:numel(subjects)
    subjectstrings{i} = sprintf('S%d',i);
end
figure(255); clf;
set(gcf,'Position',[3   912   1187 1264]);
for i=1:numel(stats_cell)
    stats = stats_cell{i};
    
    % Plot 1: Precision
    subplot(numel(stats_cell),3,(i-1)*3+1); cla; hold on;
    PlotUniqueLines([1 2], cat(1,stats(order).precision)','.',linewidth,markersize,nSubjects);
    plot(chance.precision,'k.:','linewidth',linewidth,'markersize',markersize);
    plot(best.precision,'k.--','linewidth',linewidth,'markersize',markersize);
    set(gca,'xtick',[1 2],'xticklabel',{'hBCI','hBCI+CV'},'box','on','fontname',fontname,'fontsize',fontsize)
    xlim([0.5 2.5]);
    ylim([0 100]);
    ylabel('% Precision');
    if i==numel(stats_cell)
        xlabel('Target Prediction System');
    end
    % title('Precision of Predicted Targets');
    if i==1
        legend([subjectstrings, {'Chance', 'Perfect'}],'Location','SouthEast')
    end

    % Plot 2: Recall
    subplot(numel(stats_cell),3,(i-1)*3+2); cla; hold on;
    PlotUniqueLines([1,2],cat(1,stats(order).pctFound)','.',linewidth,markersize,nSubjects);
    plot(chance.pctFound,'k.:','linewidth',linewidth,'markersize',markersize);
    plot(best.pctFound,'k.--','linewidth',linewidth,'markersize',markersize);
    % set(gca,'xtick',[1 2],'xticklabel',{'Training','TAG'},'box','on','fontname',fontname,'fontsize',fontsize)
    set(gca,'xtick',[1 2],'xticklabel',{'hBCI','hBCI+CV'},'box','on','fontname',fontname,'fontsize',fontsize)
    xlim([0.5 2.5]);
    ylim([0 100]);
    if i==numel(stats_cell)
        xlabel('Target Prediction System');
    end
    ylabel('% Targets Identified')
    title(sprintf('%s classifier',stats_names{i}));    
    % legend([subjectstrings, {'Chance'}],'Location','SouthEast')

    % Plot 3: Efficiency
    subplot(numel(stats_cell),3,(i-1)*3+3); cla; hold on;
    pctFound = cat(1,stats(order).pctFound);
    PlotUniqueLines([nan(numel(stats),1), cat(1,stats(order).pctDistance)]', [zeros(numel(stats),1), pctFound(:,2)]','.',linewidth,markersize,nSubjects);
    plot([0 chance.pctDistance],[0 100], 'k:','linewidth',linewidth,'markersize',markersize);
    plot([0 best.pctDistance],[0 100], 'k--','linewidth',linewidth,'markersize',markersize);
    set(gca,'box','on','fontname',fontname,'fontsize',fontsize)
    xlim([0 100]);
    ylim([0 100]);
    ylabel('% Targets Viewed')
    if i==numel(stats_cell)
        xlabel('% Distance Traveled');
    end    
    % title('Search Efficiency with Full System');
    % legend([subjectstrings, {'View All'}],'Location','SouthEast')
end


%% Eye Position Figures

iSubj = 9;
sess = [10 11];

EyePositionHeatMap(subjects(iSubj),sess,'3DS');
% Annotate plot
colormap gray
colormap(flipud(colormap))

iSubj = 6;
sess = [14 15]; % sessions with weird

EyePositionHeatMap(subjects(iSubj),sess,'3DS');
% Annotate plot
colormap gray
colormap(flipud(colormap))
%%
for i=1:2
    subplot(2,2,i);    
    set(gca,'xtick',[],'ytick',[],'box','on','fontname',fontname,'fontsize',fontsize);
    title(sprintf('Subject %d, Session %d',iSubj,sess(i)));
end
legend('Median eye position during session')

% PlotEyeErps_MultiSession(subjects(iSubj),sessions_cell{iSubj},0);

%% Component Figures

iSubj = 6;
EEG_noEog = RemoveEogComponents('3DS',subjects(iSubj),sessions_cell{iSubj},'-filtered-noduds',offsets(iSubj));
set(gcf,'Position',[1091        1308         284         118]);

EEG = pop_loadset('filename',sprintf('3DS-%d-all-filtered-noduds-noeog-epoched-ica.set',subjects(iSubj)));
pop_topoplot(EEG,0, [1:20] ,EEG.setname,0 ,0,'electrodes','off'); % plot scalp maps



%% FIGURES FOR RESPONSE TO REFEREES
% Learning Effects Figure
R = R_eegpsdt;
figure(265); clf;
for i=1:nSubjects
    subplot(4,3,i); cla; hold on;
    isTarg = (R(order(i)).truth>0);
    plot(find(isTarg),R(order(i)).y(isTarg),'r.');
    plot(find(~isTarg),R(order(i)).y(~isTarg),'b.');
    xlabel('trial')
    ylabel('y value')
    title(sprintf('Subject %d',i))
end
legend('target','distractor');

% Do linear fits? Anovas?

%% Make table for figure 6 (stats) results

R = R_eegpsdt;
stats = stats_eegpsdt;
tabledata = [cat(1,stats(order).precision),...
cat(1,stats(order).pctFound),...
cat(1,stats_best(order).allTourDist),...
cat(1,stats(order).pctDistance),...
cat(1,stats_best(order).bestTourDist)];

% cat(1,stats_best(order).pctDistance)];

tabledata(:,6) = tabledata(:,6)/100.*cat(1,stats_best(order).allTourDist);
tabledata(:,5:7) = tabledata(:,5:7)/1000; % put in km

% Print table headings
fprintf('Precision (%%)		Targets Identified (%%)		Distance Travelled (%%) to visit:\n');
fprintf('Subject	hBCI	hBCI + CV	hBCI	hBCI + CV	Predicted targets	True targets\n');

% Print table data
for i=1:nSubjects
%     fprintf('%d | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f\n',i,tabledata(i,:));
    fprintf('%d | %.1f | %.1f | %.1f | %.1f | %.3g | %.3g | %.3g\n',i,tabledata(i,:));
end

% fprintf('MEAN | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f\n',mean(tabledata,1));
% fprintf('STD DEV | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f\n',std(tabledata,0,1));

fprintf('MEAN | %.1f | %.1f | %.1f | %.1f | %.3g | %.3g | %.3g\n',mean(tabledata,1));
fprintf('STD DEV | %.1f | %.1f | %.1f | %.1f | %.3g | %.3g | %.3g\n',std(tabledata,0,1));

