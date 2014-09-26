% CheckPerceptualUpdate.m
%
% Look for differences in the response if the previous event was a target
% or distractor.
%
% Created 7/25/14 by DJ.

% experiment = 'sq';
% sqnumcode = 'SqNum';
% subjects = [9:11, 13:15, 17:19, 20:27];
% events = {'aD_{0/2}' 'aT_{0/2}' 'aD_{1/2}' 'aT^*_{1/2}' 'aD^*_{-/2}' 'aT_{+/2}' 'aD_{+/2}'};
% versionSuffix = 'v3pt0';

experiment = 'sf';
sqnumcode = sprintf('%s-SqNum',experiment);
subjects = [1:10 12:13];
events = {'pD_{0/2}' 'pT_{0/2}' 'pD_{1/2}' 'pT^*_{1/2}' 'pD^*_{-/2}' 'pT_{+/2}' 'pD_{+/2}'};
versionSuffix = 'v3pt2';

% experiment = 'sf3';
% sqnumcode = sprintf('%s-SqNum',experiment);
% subjects = 1:12;
% events = {'pD_{0/3}' 'pT_{0/3}' 'pD_{1/3}' 'pT_{1/3}' 'pD_{2/3}' 'pT^*_{2/3}' 'pD^*_{-/3}' 'pT_{+/3}' 'pD_{+/3}' 'pD_{-/3}'};
% versionSuffix = 'v3pt2';

types = {'Preceded by Distractor','Preceded by Target'};
nTypes = numel(types); % target or distractor
erps = nan([EEG_event.nbchan,EEG_event.pnts,nTypes,numel(subjects),numel(events)]);
nExamples = nan(nTypes,numel(subjects),numel(events));

for i = 1:numel(subjects)
    fprintf('---Subject %d/%d...\n',i,numel(subjects));
    subject = subjects(i);
    R = load(sprintf('%s-%d-GLMresults-Events-SqNum-Type-%s',experiment,subject,versionSuffix));
    
    for iEvent = 1:numel(events)
        fprintf('Event %d/%d...\n',iEvent,numel(events));
        % Get epochs
        event = events{iEvent};
        tWin = [-1 1]; % in seconds
        EEG_event = pop_epoch(R.EEG,{event},tWin);

        isPrecededByTarget = nan(1,EEG_event.trials);
        for j=1:EEG_event.trials
            lastTaskEvent = find([EEG_event.epoch(j).eventlatency{:}]<0 & ismember(EEG_event.epoch(j).eventtype,events),1,'last');     
            if ~isempty(lastTaskEvent)                
                isPrecededByTarget(j) = ~isempty(strfind(EEG_event.epoch(j).eventtype{lastTaskEvent},'T'));
            end
        end
        % Get erps
    %     isFirstHalf = SqNum<3;
    %     isSecondHalf = SqNum>3;
    %     erps(:,:,1,i) = mean(EEG_event.data(:,:,isFirstHalf),3);
    %     erps(:,:,2,i) = mean(EEG_event.data(:,:,isSecondHalf),3);    
        for k=0:1
            if any(isPrecededByTarget==k)
                erps(:,:,k+1,i,iEvent) = mean(EEG_event.data(:,:,isPrecededByTarget==k),3);
                nExamples(k+1,i,iEvent) = sum(isPrecededByTarget==k);
            end
        end
    end
end

%% Plot ERPs for same Type across SqNums
for iEvent = 1:numel(events);
    figure(iEvent); clf;
    % gridOfChans = {'F3','FZ','F4';'C3','CZ','C4';'P3','PZ','P4';'O1','OZ','O2'};
    gridOfChans = {'FZ';'CZ';'PZ';'OZ'};
%     colors = GetSquaresEventColormap({'SqNum1','SqNum2','SqNum3','SqNum4','SqNum5'});
    colors = {'b','r'};
    % PlotResponseFnsGrid(erps,{'first','second'},EEG_event.times,EEG_event.chanlocs,gridOfChans);
    PlotResponseFnsGrid(erps(:,:,:,:,iEvent),types,EEG_event.times,EEG_event.chanlocs,gridOfChans,colors);
    MakeFigureTitle([events{iEvent} ' response']);
    set(iEvent,'Position',[4 960 400 1000]);
end

CascadeFigures(1:numel(events),numel(events));
set(GetSubplots(1:numel(events)),'ylim',[-3 4],'xlim',[-500 1000]);

%% Plot ERPs for same SqNums across Types
for iType = 1:nTypes
    figure(numel(events)+iType); clf;
    % gridOfChans = {'F3','FZ','F4';'C3','CZ','C4';'P3','PZ','P4';'O1','OZ','O2'};
    gridOfChans = {'FZ';'CZ';'PZ';'OZ'};
    colors = GetSquaresEventColormap(events);
    % PlotResponseFnsGrid(erps,{'first','second'},EEG_event.times,EEG_event.chanlocs,gridOfChans);
%     PlotResponseFnsGrid(permute(erps(:,:,iSqNum,:,:),[1 2 5 4 3]),events,EEG_event.times,EEG_event.chanlocs,gridOfChans,colors);
    PlotResponseFnsGrid(permute(erps(:,:,iType,:,1:6),[1 2 5 4 3]),events(1:6),EEG_event.times,EEG_event.chanlocs,gridOfChans,colors);
    MakeFigureTitle(sprintf('Preceded-by-%s response',types{iType}));
    set(numel(events)+iType,'Position',[4 960 560 1000]);
end

CascadeFigures(numel(events)+(1:nTypes),nTypes);
set(GetSubplots(numel(events)+(1:nTypes)),'ylim',[-3 4],'xlim',[-500 1000]);

%% Plot number of examples
figure(20); cla;
% imagesc(squeeze(mean(nExamples,2))');
% set(gca,'xtick',1:nTypes,'ytick',1:numel(events),'yticklabel',events)
% colorbar
bar(squeeze(mean(nExamples,2)));
set(gca,'xtick',1:nTypes,'xticklabel',types);
xlabel('Preceding Stimulus')
ylabel('Mean nExamples')
grid on
legend(events)
set(gcf,'Position',[13 152 576 344]);

%% CALCULATE STATS
nSubjects = numel(subjects);
nT = size(erps,2);
nChans = numel(gridOfChans);
difference = nan(nChans,nT,nTypes,nTypes); % sign of significant differences between classes
disp('Getting p values...');
hwait = waitbar(0,'Getting p values...');
for l=1:numel(events)
    fprintf('Event %d/%d...\n',l,numel(events));
    for j=1:nT %nTimepoints
        waitbar(j/nT,hwait);
        for i=1:nChans   
            X = nan(nSubjects,nTypes);
            iChan = find(strcmp(gridOfChans{i},{EEG_event.chanlocs.labels}));
            for k=1:nTypes
                X(:,k) = squeeze(erps(iChan,j,l,:,k));%timecourse{i,k}(:,j);
            end
            [~,~,stats] = anova1(X,types,'off');
            comparison = multcompare(stats,'alpha',0.05,'display','off','ctype','tukey-kramer');
            sig_pos = find(comparison(:,3)>0);
            sig_neg = find(comparison(:,5)<0);
            for k=1:numel(sig_pos)
                difference(i,j,comparison(sig_pos(k),1),comparison(sig_pos(k),2)) = -1;
                difference(i,j,comparison(sig_pos(k),2),comparison(sig_pos(k),1)) = 1;
            end
            for k=1:numel(sig_neg)
                difference(i,j,comparison(sig_neg(k),1),comparison(sig_neg(k),2)) = 1;
                difference(i,j,comparison(sig_neg(k),2),comparison(sig_neg(k),1)) = -1;
            end    
        end
    end
end
close(hwait);
disp('Done!')
