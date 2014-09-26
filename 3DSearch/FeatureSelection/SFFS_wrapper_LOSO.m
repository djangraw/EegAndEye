function sffs = SFFS_wrapper_LOSO(subjects,sessions_cell,offsets)

% Run TEMP_ClassifyWithPupilSize_v2 first, with useEEG turned on. this will
% get you variables subjects,sessions_cell,offsets, R(i).y_EEG,ps_all,t_epoch,binstart,binwidth
% 
% Created 10/2/13 by DJ based on SFFS_wrapper.m.

%% Get features
[dwell_time, sac_size, sac_speed, fixtime_max, fixtime_fixMean, nsac, fixtime_first, rxn_time, fixtime_b4, fixtime_after, sacsize_away, dist2obj_b4, dist2obj_first, dist2obj_last, dist2obj_after, fixtime_last] = GetVariousEyeFeatures(subjects,sessions_cell);
% bigEyeFeature = AppendFeatures(dwell_time, sac_size, sac_speed, fixtime_max, fixtime_fixMean, nsac, fixtime_first);
bigEyeFeature = AppendFeatures(dwell_time, sac_size, fixtime_mean, nsac, fixtime_first, rxn_time, fixtime_b4, fixtime_after, sacsize_away, dist2obj_b4, dist2obj_first, dist2obj_last, dist2obj_after, fixtime_last);

%% Get pupil features
[ps_bin, ps_max, ps_latency, ps_deriv] = GetVariousPupilFeatures(ps_all,t_epoch,binstart,binwidth);
bigPupFeature = AppendFeatures(ps_bin, ps_max, ps_latency, ps_deriv);

%% Combine features
bigFeature = AppendFeatures(bigEyeFeature,bigPupFeature);
for iSubj = 1:numel(bigFeature)   
    bigFeature{iSubj} = cat(2,R(iSubj).y_EEG,bigFeature{iSubj});
end
nFeats = size(bigFeature{1},2);

%% Set up
[EEGDATA,YDATA] = deal(cell(1,numel(subjects)));
disp('Loading Data...')
for iSubj = 1:numel(subjects)
    load(sprintf('ALLEEG%d_NoeogEpochedIcaCropped.mat',subjects(iSubj))); % ALLEEG
    y = loadBehaviorData(subjects(iSubj),sessions_cell{iSubj},'3DS');
    EEGDATA{iSubj} = ALLEEG;
    YDATA{iSubj} = y;
end
disp('Done!')


%% SET UP SFFS ACROSS SUBJECTS
nSubjects = numel(subjects);
NumFeatComb = nFeats;
CostFunction = 'ClassifierWrapper_forSFFS_xsubj';
cvmode = '10fold';

% Set up features
otherInputs = {EEGDATA, YDATA, offsets};
featcell = cell(nFeats,nSubjects);
for i=1:nFeats    
    for iSubj = 1:nSubjects
        featcell{i,iSubj} = bigFeature{iSubj}(:,i);
    end
end
        
% Get truth values
truth = cell(1,nSubjects);
for iSubj = 1:nSubjects
    truth{iSubj} = [R(iSubj).truth]';
end

% % find cross-validation trials
% nFolds = 10;
% [incTrials,valTrials] = deal(cell(nSubjects,nFolds));
% [truth_train, truth_test] = deal(cell(size(truth)));
% for iSubj = 1:nSubjects
%     nTrials = length(featcell{1,iSubj});
%     cv = setCrossValidationStruct(cvmode,nTrials);
%     incTrials(iSubj,:) = cv.incTrials;
%     valTrials(iSubj,:) = cv.valTrials;
%     for iFold=1:nFolds
%         truth_train{iSubj,iFold} = truth{iSubj}(incTrials{iSubj,iFold});
%         truth_test{iSubj,iFold} = truth{iSubj}(valTrials{iSubj,iFold});
%     end
% end

%% RUN TRAINING: leave-one-subj-out xval (score = avg. 10-fold AUC across subjects,
%% leave out all of 1 subject's data)
[featsBest, scoreBest, featsAll, scoreAll] = deal(cell(1,nSubjects));
fprintf('Started xval training at time: %s \n',datestr(now));
nFolds = nSubjects;
parfor iFold = 1:nSubjects
    if ~isempty(featsBest{iFold})
        fprintf('Skipping fold %d...\n',iFold);
        continue
    end
    tstart = tic;
    fprintf('Starting fold %d...\n',iFold);
    featcell_cropped = featcell(:,setdiff(1:nSubjects,iFold)); % leave out subj iFold
    truth_cropped = truth(:,setdiff(1:nSubjects,iFold)); % leave out subj iFold
    [featsBest{iFold}, scoreBest{iFold}, featsAll{iFold}, scoreAll{iFold}] = SequentialForwardFloatingSelection_DJ(featcell_cropped,featcell_cropped,CostFunction,NumFeatComb,{truth_cropped});
%     [featsBest{iFold}, scoreBest{iFold}, featsAll{iFold}, scoreAll{iFold}] = SequentialForwardFloatingSelection_DJ_ADDON(featcell_cropped,featcell_cropped,CostFunction,NumFeatComb,featsAll{iFold}, scoreAll{iFold}, {truth_train(:,iFold)});
    % DO TESTING LATER.
    fprintf('Fold %d done! took %g seconds.\n',iFold,toc(tstart));
end
fprintf('Done with xval training at time: %s \n',datestr(now));


%% RUN TESTING: 10-fold xval on LEFT-OUT TRIALS (10% of original set)
[scoreTest] = deal(cell(1,nFolds));
fprintf('Started xval testing at time: %s \n',datestr(now));
for iFold = 1:nFolds
    if ~isempty(scoreTest{iFold})
        fprintf('Skipping fold %d...\n',iFold);
        continue
    end
    tstart = tic;
    fprintf('Starting fold %d...\n',iFold);
    featcell_test = featcell(:,iFold); % use only subj iFold
    truth_test = truth(:,iFold);

    scoreTest{iFold} = nan(1,NumFeatComb);
    for j=1:NumFeatComb % # features to include in set
        if ~isempty(sffs_xsubj.features{iFold}{j})
            featcell_test_cropped = featcell_test(sffs_xsubj.features{iFold}{j},:);
            scoreTest{iFold}(j) = ClassifierWrapper_forSFFS_xsubj(featcell_test_cropped,featcell_test_cropped,truth_test);
        else
            scoreTest{iFold}(j) = -inf;
        end
    end
        
    fprintf('Fold %d done! took %g seconds.\n',iFold,toc(tstart));
end
fprintf('Done with xval testing at time: %s \n',datestr(now));


%% Compile results and SAVE
sffs_xsubj.features = featsAll;
sffs_xsubj.scores = scoreAll;

% Get feature names and orders
featurenames = {'EEG(150)', 'EEG(250)', 'EEG(350)', 'EEG(450)', 'EEG(550)', 'EEG(650)', 'EEG(750)', 'EEG(850)', 'EEG(950)', ...
    'dwellTime','size_sacTo','dur_fixMean','nFix', 'dur_fixFirst','latency_fixFirst', 'dur_fixBefore', 'dur_fixAfter', 'size_sacAway','dist_fixBefore', 'dist_fixFirst', 'dist_fixLast', 'dist_fixAfter','dur_fixLast'...
    'PD(250)','PD(750)','PD(1250)','PD(1750)','PD(2250)','PD(2750)','PDmax','latency_PDmax','PD''(250)','PD''(750)','PD''(1250)','PD''(1750)','PD''(2250)','PD''(2750)'};

NEWfeaturenames = {'EEG(150)', 'EEG(250)', 'EEG(350)', 'EEG(450)', 'EEG(550)', 'EEG(650)', 'EEG(750)', 'EEG(850)', 'EEG(950)', ...
    'nFix', 'dist_fixBefore','size_sacTo','dist_fixFirst', 'dist_fixLast','size_sacAway','dist_fixAfter', 'dur_fixBefore','latency_fixFirst','dur_fixFirst','dur_fixMean','dwellTime','dur_fixLast','dur_fixAfter',...
    'latency_PDmax','PDmax','PD(250)','PD(750)','PD(1250)','PD(1750)','PD(2250)','PD(2750)','PD''(250)','PD''(750)','PD''(1250)','PD''(1750)','PD''(2250)','PD''(2750)'};
for i=1:nFeats
    featorder(i) = find(strcmp(NEWfeaturenames{i},featurenames));
end
sffs_xsubj.featurenames = featurenames;
sffs_xsubj.featorder = featorder;
sffs_xsubj.scoreTest = cat(1,scoreTest{:});

% Save result
fprintf('Saving results as TEMP_SFFS_loso_%dfeats.mat...\n',nFeats);
save(sprintf('TEMP_SFFS_loso_%dfeats.mat',nFeats),'sffs_xsubj');
disp('Done!')

%% Get features of size X

% Set plot parameters
fontname = 'Helvetica'; % Futura, Helvetica, Gotham, Georgia are all good
fontsize = 12;
linewidth = 2;
markersize = 20;


iEEG = 1:9;
iNum = 10;
iDist = 11:16;
iTime = 17:24;
iPD = 25:31;
idPDdt = 32:37;

nFeats = size(bigFeature{1},2);
isInSet = zeros(NumFeatComb,nFeats,numel(subjects));
for iSubj=1:numel(subjects)
    thisset = [];
    for i=1:NumFeatComb
%         if ~isempty(sffs.features{iSubj}{i});
            thisset = sffs_xsubj.features{iSubj}{i};
%         end
        isInSet(i,thisset,iSubj) = 1;
    end
end
figure(333); clf;
set(gcf,'Position',[4   953   996   531]);
set(gca,'box','on','fontname',fontname,'fontsize',fontsize);
imagesc(mean(isInSet(:,featorder,:),3));
colorbar('fontname',fontname,'fontsize',fontsize);
colormap gray
% colormap(flipud(colormap('gray')))
set(gca,'xtick',1:nFeats,'xticklabel',show_symbols(NEWfeaturenames),'ydir','normal')
set(gca,'Position',[0.100    0.2300    0.800    0.7050]);
% foo = rotateticklabel(gca,90);
% set(foo,'fontname',fontname,'fontsize',fontsize);
ylabel('Set size')
title('SFFS use rate as a function of set size')

% Plot usage at certain # of features
iFeatComb = 9;%11;%NumFeatComb; % pick # of features to use
figure(334); clf; hold on
set(gcf,'Position',[8   343   991   539]);   
set(gca,'box','on','fontname',fontname,'fontsize',fontsize);
useRate = mean(isInSet(iFeatComb,featorder,:),3);
bar(iEEG,useRate(iEEG),'b');
bar(iTime,useRate(iTime),'r');
bar(iDist,useRate(iDist),'facecolor',[0 .5 0]);
bar(iPD,useRate(iPD),'facecolor',[.5 0 1]);
bar(iNum,useRate(iNum),'c');
bar(idPDdt,useRate(idPDdt),'facecolor',[1 .5 0]);
% bar(useRate);
set(gca,'xtick',1:nFeats,'xticklabel',show_symbols(NEWfeaturenames),'ydir','normal')
set(gca,'Position',[0.100    0.2300    0.850    0.7050]);
xlim([0 nFeats+1])
foo = rotateticklabel(gca,90);
set(foo,'fontname',fontname,'fontsize',fontsize);
ylabel('use rate')
title(sprintf('Use rate at set size = %d',iFeatComb))

%% Plot Az values
allscores = cat(1,sffs_xsubj.scores{:});
meanscore = mean(allscores(:,1:NumFeatComb),1);
stderrscore = std(allscores(:,1:NumFeatComb),[],1)/sqrt(nFolds);
figure(335); clf; hold on
set(gca,'box','on','fontname',fontname,'fontsize',fontsize);
ErrorPatch(0:NumFeatComb,[0.5 meanscore],[0 stderrscore],'b','b');

% Add testing Az's
allscores = sffs_xsubj.scoreTest;
meanscore = mean(allscores(:,1:NumFeatComb),1);
stderrscore = std(allscores(:,1:NumFeatComb),[],1)/sqrt(nFolds);
% Plot
ErrorPatch(0:NumFeatComb,[0.5 meanscore],[0 stderrscore],[0 .5 0],[0 .5 0]);

% Annotate plot
xlabel('Set size');
ylabel('AUC')
plot([0 nFeats+1],[0.5 0.5],'k--','LineWidth',linewidth)
ylim([0.3 1]);
xlim([0 nFeats]);

%% extra: add non-LOSO testing
figure(335); clf; hold on
set(gca,'box','on','fontname',fontname,'fontsize',fontsize);
foo = load(sprintf('TEMP_SFFS_xsubj_%dfeats.mat',nFeats));

% Plot LOSO training
allscores = cat(1,sffs_xsubj.scores{:});
meanscore = mean(allscores(:,1:NumFeatComb),1);
stderrscore = std(allscores(:,1:NumFeatComb),[],1)/sqrt(nFolds);
ErrorPatch(0:NumFeatComb,[0.5 meanscore],[0 stderrscore],'g','g');

% Plot LOSO testing Az's
allscores = sffs_xsubj.scoreTest;
meanscore = mean(allscores(:,1:NumFeatComb),1);
stderrscore = std(allscores(:,1:NumFeatComb),[],1)/sqrt(nFolds);
ErrorPatch(0:NumFeatComb,[0.5 meanscore],[0 stderrscore],[0 .5 0],[0 .5 0]);

% Plot training
allscores = cat(1,foo.sffs_xsubj.scores{:});
meanscore = mean(allscores(:,1:NumFeatComb),1);
stderrscore = std(allscores(:,1:NumFeatComb),[],1)/sqrt(nFolds);
ErrorPatch(0:NumFeatComb,[0.5 meanscore],[0 stderrscore],'b','b');

% Plot testing
allscores = foo.sffs_xsubj.scoreTest;
meanscore = mean(allscores(:,1:NumFeatComb),1);
stderrscore = std(allscores(:,1:NumFeatComb),[],1)/sqrt(nFolds);
ErrorPatch(0:NumFeatComb,[0.5 meanscore],[0 stderrscore],'r','r');

% Annotate plot
xlabel('Set size');
ylabel('AUC')
plot([0 nFeats+1],[0.5 0.5],'k--','LineWidth',linewidth)
ylim([0.3 1]);
xlim([0 nFeats]);

MakeLegend([0 0 1; 1 0 0; 0 1 0; 0 0.5 0],{'10-fold Training','10-fold Testing','LOSO Training','LOSO Testing'},2);


%% Plot set size at which usage > 0.5
figure(332); clf; hold on
set(gcf,'Position',[1000   343   991   539]);   
set(gca,'box','on','fontname',fontname,'fontsize',fontsize);
useRate = mean(isInSet(:,featorder,:),3);
useRate = cat(1,zeros(1,size(useRate,2)),useRate);
% for i=1:nFeats
%     over50size(i) = find(useRate(:,i)>0.5,1);
% end
% [~,over50order] = sort(over50size,'ascend');
% orderedsize = over50size(over50order);

% bar(find(ismember(over50order,iEEG)),orderedsize(ismember(over50order,iEEG)),'b');
% bar(find(ismember(over50order,iTime)),orderedsize(ismember(over50order,iTime)),'r');
% bar(find(ismember(over50order,iDist)),orderedsize(ismember(over50order,iDist)),'facecolor',[0 .5 0]);
% bar(find(ismember(over50order,iPD)),orderedsize(ismember(over50order,iPD)),'facecolor',[.5 0 1]);
% bar(find(ismember(over50order,iNum)),orderedsize(ismember(over50order,iNum)),'c');
% bar(find(ismember(over50order,idPDdt)),orderedsize(ismember(over50order,idPDdt)),'facecolor',[1 .5 0]);

bar(iEEG,mean(useRate(:,iEEG)),'b');
bar(iTime,mean(useRate(:,iTime)),'r');
bar(iDist,mean(useRate(:,iDist)),'facecolor',[0 .5 0]);
bar(iPD,mean(useRate(:,iPD)),'facecolor',[.5 0 1]);
bar(iNum,mean(useRate(:,iNum)),'c');
bar(idPDdt,mean(useRate(:,idPDdt)),'facecolor',[1 .5 0]);
chanceUse = mean(useRate(:));
plot([0 nFeats+1],[chanceUse, chanceUse],'k--','linewidth',2);

% bar(iEEG,over50size(iEEG),'b');
% bar(iTime,over50size(iTime),'r');
% bar(iDist,over50size(iDist),'facecolor',[0 .5 0]);
% bar(iPD,over50size(iPD),'facecolor',[.5 0 1]);
% bar(iNum,over50size(iNum),'c');
% bar(idPDdt,over50size(idPDdt),'facecolor',[1 .5 0]);

% set(gca,'xtick',1:nFeats,'xticklabel',show_symbols(NEWfeaturenames(over50order)),'ydir','normal')
set(gca,'xtick',1:nFeats,'xticklabel',show_symbols(NEWfeaturenames),'ydir','normal')
set(gca,'Position',[0.100    0.2300    0.850    0.7050]);
xlim([0 nFeats+1])
foo = rotateticklabel(gca,90);
set(foo,'fontname',fontname,'fontsize',fontsize);
% ylabel('Set size at which use rate > 0.5')
ylabel('Mean use rate across all set sizes')
% title('Feature ranking')
