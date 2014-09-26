function sffs = SFFS_wrapper(subjects,sessions_cell,offsets)


%% Get features
[dwell_time, sac_size, sac_speed, fixtime_max, fixtime_mean, nsac, fixtime_first, rxn_time, fixtime_b4, fixtime_after, sacsize_away, dist2obj_b4, dist2obj_first, dist2obj_last, dist2obj_after, fixtime_last] = GetVariousEyeFeatures(subjects,sessions_cell);
% bigEyeFeature = AppendFeatures(dwell_time, sac_size, sac_speed, fixtime_max, fixtime_mean, nsac, fixtime_first);
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
%% RUN SFFS FOR EACH SUBJECT
NumFeatComb = size(bigFeature{1},2);
CostFunction = 'ClassifierWrapper_forSFFS';
[featsBest, scoreBest, featsAll, scoreAll] = deal(cell(1,numel(subjects)));
parfor iSubj = 1:numel(subjects)
    try
        % Run SFFS
        fprintf('====== Subject %d ======\n',subjects(iSubj));
    %     otherInputs = {ALLEEG, y, offsets(iSubj)};
        otherInputs = {EEGDATA{iSubj}, YDATA{iSubj}, offsets(iSubj)};
        [featsBest{iSubj}, scoreBest{iSubj}, featsAll{iSubj}, scoreAll{iSubj}] = SequentialForwardFloatingSelection_DJ(bigFeature{iSubj}',bigFeature{iSubj}',CostFunction,NumFeatComb,otherInputs);
    catch
        msgbox(sprintf('SFFS Failed for Subject %d!',subjects(iSubj)));
    end
end

% Compile results
sffs.features = featsAll;
sffs.scores = scoreAll;

%% Get features of size X

isInSet = zeros(NumFeatComb,nFeats,numel(subjects));
for iSubj=1:numel(subjects)
    thisset = [];
    for i=1:NumFeatComb
%         if ~isempty(sffs.features{iSubj}{i});
            thisset = sffs.features{iSubj}{i};
%         end
        isInSet(i,thisset,iSubj) = 1;
    end
end
figure(333); clf;
imagesc(mean(isInSet,3));
colorbar;
featurenames = {'EEG(100)', 'EEG(200)', 'EEG(300)', 'EEG(400)', 'EEG(500)', 'EEG(600)', 'EEG(700)', 'EEG(800)', 'EEG(900)', ...
    'dt','sacsize','mean tfix','nfix', 'first tfix','RT', 'fixtime b4', 'fixtime after', 'sacsize away','dist2obj before', 'dist2obj first', 'dist2obj last'...
    'PD(0)','PD(500)','PD(1000)','PD(1500)','PD(2000)','PD(2500)','PD max','PD latency','PD''(0)','PD''(500)','PD''(1000)','PD''(1500)','PD''(2000)','PD''(2500)'};


set(gca,'xtick',1:nFeats,'xticklabel',featurenames,'ydir','normal')
ylabel('Set size')

figure(334); clf;
iFeatComb = NumFeatComb;
useRate = mean(isInSet(iFeatComb,:,:),3);
plot(useRate,'.-');
set(gca,'xtick',1:nFeats,'xticklabel',featurenames,'ydir','normal')
ylabel('use rate')




%% SET UP SFFS ACROSS SUBJECTS
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

% find cross-validation trials
nFolds = 10;
[incTrials,valTrials] = deal(cell(nSubjects,nFolds));
[truth_train, truth_test] = deal(cell(size(truth)));
for iSubj = 1:nSubjects
    nTrials = length(featcell{1,iSubj});
    cv = setCrossValidationStruct(cvmode,nTrials);
    incTrials(iSubj,:) = cv.incTrials;
    valTrials(iSubj,:) = cv.valTrials;
    for iFold=1:nFolds
        truth_train{iSubj,iFold} = truth{iSubj}(incTrials{iSubj,iFold});
        truth_test{iSubj,iFold} = truth{iSubj}(valTrials{iSubj,iFold});
    end
end

%% RUN TRAINING: 10-fold xval (score = avg. 10-fold AUC across subjects, leave out 10% of each subject's data)
[featsBest, scoreBest, featsAll, scoreAll] = deal(cell(1,numel(subjects)));
fprintf('Started 10-fold at time: %s \n',datestr(now));
parfor iFold = 1:nFolds
    if ~isempty(featsBest{iFold})
        fprintf('Skipping fold %d...\n',iFold);
        continue
    end
    tstart = tic;
    fprintf('Starting fold %d...\n',iFold);
    featcell_cropped = cell(size(featcell));    
    for i=1:nFeats
        for iSubj = 1:nSubjects
            featcell_cropped{i,iSubj} = featcell{i,iSubj}(incTrials{iSubj,iFold});            
        end
    end    
    [featsBest{iFold}, scoreBest{iFold}, featsAll{iFold}, scoreAll{iFold}] = SequentialForwardFloatingSelection_DJ(featcell_cropped,featcell_cropped,CostFunction,NumFeatComb,{truth_train(:,iFold)});
%     [featsBest{iFold}, scoreBest{iFold}, featsAll{iFold}, scoreAll{iFold}] = SequentialForwardFloatingSelection_DJ_ADDON(featcell_cropped,featcell_cropped,CostFunction,NumFeatComb,featsAll{iFold}, scoreAll{iFold}, {truth_train(:,iFold)});
    % DO TESTING LATER.
    fprintf('Fold %d done! took %g seconds.\n',iFold,toc(tstart));
end
fprintf('Done with 10-fold at time: %s \n',datestr(now));

%% Compile results
sffs_xsubj.features = featsAll;
sffs_xsubj.scores = scoreAll;

%% Get features of size X

% Set plot parameters
fontname = 'Helvetica'; % Futura, Helvetica, Gotham, Georgia are all good
fontsize = 12;
linewidth = 2;
markersize = 20;
% Get feature names and orders
featurenames = {'EEG(150)', 'EEG(250)', 'EEG(350)', 'EEG(450)', 'EEG(550)', 'EEG(650)', 'EEG(750)', 'EEG(850)', 'EEG(950)', ...
    'dwellTime','sizeSac_to','durFix_mean','nFix', 'durFix_first','reactionTime', 'durFix_before', 'durFix_after', 'sizeSac_away','dist_before', 'dist_first', 'dist_last', 'dist_after','durFix_last'...
    'PD(250)','PD(750)','PD(1250)','PD(1750)','PD(2250)','PD(2750)','PDmax','latency_PDmax','PD''(250)','PD''(750)','PD''(1250)','PD''(1750)','PD''(2250)','PD''(2750)'};

NEWfeaturenames = {'EEG(150)', 'EEG(250)', 'EEG(350)', 'EEG(450)', 'EEG(550)', 'EEG(650)', 'EEG(750)', 'EEG(850)', 'EEG(950)', ...
    'nFix', 'dist_before','sizeSac_to','dist_first', 'dist_last','sizeSac_away','dist_after', 'reactionTime','dwellTime','durFix_before','durFix_first','durFix_mean','durFix_last','durFix_after',...
    'latency_PDmax','PDmax','PD(250)','PD(750)','PD(1250)','PD(1750)','PD(2250)','PD(2750)','PD''(250)','PD''(750)','PD''(1250)','PD''(1750)','PD''(2250)','PD''(2750)'};
for i=1:nFeats
    featorder(i) = find(strcmp(NEWfeaturenames{i},featurenames));
end
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
iFeatComb = 11;%NumFeatComb; % pick # of features to use
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

% Plot Az values
allscores = cat(1,sffs_xsubj.scores{:});
meanscore = mean(allscores(:,1:NumFeatComb),1);
stderrscore = std(allscores(:,1:NumFeatComb),[],1)/sqrt(nFolds);
figure(335); clf; hold on
set(gca,'box','on','fontname',fontname,'fontsize',fontsize);
ErrorPatch(0:NumFeatComb,[0.5 meanscore],[0 stderrscore],'b','b');
% errorbar(1:NumFeatComb, meanscore,stderrscore);
xlabel('# features included');
ylabel('Training AUC')
plot([0 nFeats+1],[0.5 0.5],'k--','LineWidth',linewidth)
ylim([0.3 1]);
xlim([0 nFeats+1]);

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
