% TEMP_MakeEfficiencyPlots.m
%
% Created 8/20/12 by DJ for EMBC 2012 conference.
% Updated 3/29/13 by DJ - subjectstring for variable # of subjects


%% SET UP
subjects = 22:29;
% subjects = [18 19 20 22 23];
sessions = {2:14, [3 6:17] 1:15, 1:15, 1:15, 1:15, 1:15, 1:15};
% sessions = {1:10, 2:11, 1:12, 2:14, [3 6:17]};
folders = {'2013-03-27-Pilot-S22','2013-03-28-Pilot-S23' '2013-04-29-Pilot-S24','2013-05-01-Pilot-S25' '2013-05-02-Pilot-S26','2013-05-03-Pilot-S27' '2013-05-07-Pilot-S28', '2013-05-10-Pilot-S29'};
% folders = {'2011-10-25-Pilot', '2011-11-08-Pilot','2012-01-12-Pilot','2013-03-27-Pilot-S22','2013-03-28-Pilot-S23'};
version = 'filtered-noduds-ica';
% version = 'nofrontal-noduds';
Numbers = GetNumbers; 
eventnumbers = [Numbers.SACCADE_TO+Numbers.TARGET, Numbers.SACCADE_TO+Numbers.DISTRACTOR];
eventnames = {'targsac','distsac'};
usegridconstraints = true;
cvmode = '10fold';
useica = true;
basedir = '/Users/dave/Documents/Data/3DSearch';

%% CALCULATE
[stats,fwdModel,v,iObjects_eeg_pt,chanlocs] = deal(cell(1,3));
cd(basedir)
for i=1:numel(subjects)
    cd(folders{i});
    filename = sprintf('3DS-%d-all-%s.set',subjects(i),version); % use suffix indicated by input string 'version'
    [ALLEEG, EEG, CURRENTSET] = EpochEeglabFile(filename,[-1000 2000],[-200 0],eventnumbers,eventnames,'None');
    [iObjects_tag_pt, iObjects_eeg_pt{i},fwdModel{i},v{i}] = RunResultsThroughTAG(subjects(i),sessions{i},ALLEEG([3 2]),cvmode,useica);
    [stats{i},TspTour] = CalculateImprovement(subjects(i),sessions{i},'GridHuge.png',iObjects_eeg_pt{i},usegridconstraints);
    chanlocs{i} = ALLEEG(1).chanlocs;
    cd(basedir);
end
stats = [stats{:}];
%% SET CHANCE LEVELS

chance.precision = [25 25];
chance.pctFound = [20/56*5 25];
chance.pctDistance = 100;

%% PLOT
% Set plot parameters
fontname = 'Futura';
fontsize = 15;
linewidth = 2;
markersize = 20;

subjectstrings = cell(1,numel(subjects));
for i=1:numel(subjects)
    subjectstrings{i} = sprintf('S%d',i);
end

% Make plots
figure(251); clf;
cla; hold on;
plot(cat(1,stats.precision)', '.-','linewidth',linewidth,'markersize',markersize);
plot(chance.precision,'k.--','linewidth',linewidth,'markersize',markersize);
set(gca,'xtick',[1 2],'xticklabel',{'EEG','EEG-TAG'},'box','on','fontname',fontname,'fontsize',fontsize)
xlim([0 3]);
ylim([0 100]);
ylabel('% Precision');
xlabel('Target Prediction System');
title('Precision of Predicted Targets');
legend([subjectstrings, {'Chance'}],'Location','SouthEast')

figure(252); clf;
cla; hold on;
plot(cat(1,stats.pctFound)', '.-','linewidth',linewidth,'markersize',markersize);
plot(chance.pctFound,'k.--','linewidth',linewidth,'markersize',markersize);
set(gca,'xtick',[1 2],'xticklabel',{'Training','TAG'},'box','on','fontname',fontname,'fontsize',fontsize)
xlim([0 3]);
ylim([0 100]);
ylabel('% Targets Identified')
title('Percentage of True Targets Identified');
legend([subjectstrings, {'Chance'}],'Location','SouthEast')

figure(253); clf;
cla; hold on;
pctFound = cat(1,stats.pctFound);
plot([zeros(numel(stats),1), cat(1,stats.pctDistance)]', [zeros(numel(stats),1), pctFound(:,2)]', '.-','linewidth',linewidth,'markersize',markersize);
plot([0 100], [0 chance.pctDistance],'k.--','linewidth',linewidth,'markersize',markersize);
set(gca,'box','on','fontname',fontname,'fontsize',fontsize)
xlim([0 100]);
ylim([0 100]);
ylabel('% Targets Viewed')
xlabel('% Distance traveled (normalized)')
title('Search Efficiency with Full System');
legend([subjectstrings, {'View All'}],'Location','SouthEast')

%% Plot average forward models
bintimes = 0:100:1000;
figure(255); clf;
channies = [chanlocs{:}];
chans = unique({channies.labels});
for i=1:numel(subjects)
    chans = chans(ismember(chans,{chanlocs{i}.labels}));
end
chanlocs_all = chanlocs{1}(ismember({chanlocs{1}.labels},chans));
% chanmatrix = [1:10, 0; 0 0, 11:17, 0, 0; 18:27, 0; 28:38; 39:40,0,41:48; 0,0,0,0,49:51,0,0,0,0; 0,0,0,0,52:54,0,0,0,0];

fwds = zeros(numel(chans), size(fwdModel{1},2),numel(subjects));
twts = zeros(size(v{1},1), size(v{1},2),numel(subjects));
for i=1:numel(subjects)
   fwds(:,:,i) = mean(fwdModel{i}(ismember({chanlocs{i}.labels},chans),:,:),3);
   twts(:,:,i) = mean(v{i},3);
end

% cmax = max(max(abs(mean(fwds,3))));
cmax = 2;
for i=1:size(fwds,2)
    subplot(4,5,rem(i-1,2)*5+floor((i-1)/2)+1,'fontname',fontname,'fontsize',fontsize)
    topoplot(mean(fwds(:,i,:),3),chanlocs_all,'maplimits',[-cmax cmax],'plotrad',0.5);    
%     topoplot([],chanlocs{1}(ismember({chanlocs{1}.labels},chans)),'plotrad',0.5);
%     topoplot([],chanlocs_all,'electrodes','off','headrad',0.63);
%     topoplot(mean(fwds(:,i,:),3), chanlocs_all,'plotgrid',chanmatrix);    
%     title(sprintf('%d-%dms',bintimes(i),bintimes(i+1)));
end
axes('Position',[.85 .55 .1 .24],'CLim',[-cmax cmax],'visible','off','fontname',fontname,'fontsize',fontsize);
colorbar('fontname',fontname,'fontsize',fontsize)

h = MakeFigureTitle('Forward Models (Mean Across Subjects) (uV)');
set(h,'fontname',fontname,'fontsize',fontsize,'fontweight','normal')

subplot(2,1,2,'fontname',fontname,'fontsize',fontsize); hold on;
xbin = bintimes(1:end-1)+diff(bintimes)/2;
plot(xbin,permute(twts,[2 3 1]),'.-','linewidth',linewidth,'markersize',markersize);
legend('S1', 'S2', 'S3', 'Location','SouthWest')
% plot(xbin,mean(twts,3),'.-','linewidth',linewidth,'markersize',markersize
% );
% for i=1:numel(xbin)
%     plot([xbin(i),xbin(i)],[mean(twts(:,i,:),3)-std(twts(:,i,:),[],3)/sqrt(size(twts,3)), ...
%         mean(twts(:,i,:),3)+std(twts(:,i,:),[],3)/sqrt(size(twts,3))],...
%         '--','linewidth',linewidth,'markersize',markersize)
% end
set(gca,'xgrid','on','box','on')
ylabel('temporal weights (A.U.)')
xlabel('time between saccade offset and bin center (ms)')
