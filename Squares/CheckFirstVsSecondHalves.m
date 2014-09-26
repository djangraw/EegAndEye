% CheckFirstVsSecondHalves
%
% Created 7/22/14 by DJ.

experiment = 'sq';
sqnumcode = 'SqNum';
subjects = [9:11, 13:15, 17:19, 20:27];
events = {'aD_{0/2}' 'aT_{0/2}' 'aD_{1/2}' 'aT^*_{1/2}' 'aD^*_{-/2}' 'aT_{+/2}' 'aD_{+/2}'};
versionSuffix = 'v3pt0';

% experiment = 'sf';
% sqnumcode = sprintf('%s-SqNum',experiment);
% subjects = [1:10 12:13];
% events = {'pD_{0/2}' 'pT_{0/2}' 'pD_{1/2}' 'pT^*_{1/2}' 'pD^*_{-/2}' 'pT_{+/2}' 'pD_{+/2}'};
% versionSuffix = 'v3pt2';

% experiment = 'sf3';
% sqnumcode = sprintf('%s-SqNum',experiment);
% subjects = 1:12;
% events = {'pD_{0/3}' 'pT_{0/3}' 'pD_{1/3}' 'pT_{1/3}' 'pD_{2/3}' 'pT^*_{2/3}' 'pD^*_{-/3}' 'pT_{+/3}' 'pD_{+/3}' 'pD_{-/3}'};
% versionSuffix = 'v3pt2';

erps = nan([EEG_event.nbchan,EEG_event.pnts,5,numel(subjects),numel(events)]);
nExamples = nan(5,numel(subjects),numel(events));

for i = 1:numel(subjects)
    fprintf('Subject %d/%d...\n',i,numel(subjects));
    subject = subjects(i);
    R = load(sprintf('%s-%d-GLMresults-Events-SqNum-Type-%s',experiment,subject,versionSuffix));

    for iEvent = 1:numel(events)
        fprintf('Event %d/%d...\n',iEvent,numel(events));
        % Get epochs
        event = events{iEvent};
        tWin = [-.5 1]; % in seconds
        EEG_event = pop_epoch(R.EEG,{event},tWin);

        SqNum = nan(1,EEG_event.trials);
        for j=1:EEG_event.trials
            sqnum_anchor = find([EEG_event.epoch(j).eventlatency{:}]==0 & strncmp(sqnumcode,EEG_event.epoch(j).eventtype,length(sqnumcode)));     
            if ~isempty(sqnum_anchor)
                SqNum(j) = str2num(EEG_event.epoch(j).eventtype{sqnum_anchor}(end));
            end
        end
        % Get erps
    %     isFirstHalf = SqNum<3;
    %     isSecondHalf = SqNum>3;
    %     erps(:,:,1,i) = mean(EEG_event.data(:,:,isFirstHalf),3);
    %     erps(:,:,2,i) = mean(EEG_event.data(:,:,isSecondHalf),3);    
        for k=1:5
            if any(SqNum==k)
                erps(:,:,k,i,iEvent) = mean(EEG_event.data(:,:,SqNum==k),3);
                nExamples(k,i,iEvent) = sum(SqNum==k);
            end
        end
    end
end

%% Plot ERPs for same Type across SqNums
types = {'SqNum1','SqNum2','SqNum3','SqNum4','SqNum5'};
nTypes = numel(types);
for iEvent = 1:numel(events);
    figure(iEvent); clf;
    % gridOfChans = {'F3','FZ','F4';'C3','CZ','C4';'P3','PZ','P4';'O1','OZ','O2'};
    gridOfChans = {'FZ';'CZ';'PZ';'OZ'};
    colors = GetSquaresEventColormap(types);
    % PlotResponseFnsGrid(erps,{'first','second'},EEG_event.times,EEG_event.chanlocs,gridOfChans);
    PlotResponseFnsGrid(erps(:,:,:,:,iEvent),types,EEG_event.times,EEG_event.chanlocs,gridOfChans,colors);
    MakeFigureTitle([events{iEvent} ' response']);
    set(iEvent,'Position',[4 960 400 1000]);
end

CascadeFigures(1:numel(events),numel(events));
set(GetSubplots(1:numel(events)),'ylim',[-3 4]);

%% Plot ERPs for same SqNums across Types
for iSqNum = 1:5
    figure(numel(events)+iSqNum); clf;
    % gridOfChans = {'F3','FZ','F4';'C3','CZ','C4';'P3','PZ','P4';'O1','OZ','O2'};
    gridOfChans = {'FZ';'CZ';'PZ';'OZ'};
    colors = GetSquaresEventColormap(events);
    % PlotResponseFnsGrid(erps,{'first','second'},EEG_event.times,EEG_event.chanlocs,gridOfChans);
%     PlotResponseFnsGrid(permute(erps(:,:,iSqNum,:,:),[1 2 5 4 3]),events,EEG_event.times,EEG_event.chanlocs,gridOfChans,colors);
    PlotResponseFnsGrid(permute(erps(:,:,iSqNum,:,1:4),[1 2 5 4 3]),events(1:4),EEG_event.times,EEG_event.chanlocs,gridOfChans,colors);
    MakeFigureTitle(sprintf('SqNum%d response',iSqNum));
    set(numel(events)+iSqNum,'Position',[4 960 560 1000]);
end

CascadeFigures(numel(events)+(1:5),5);
set(GetSubplots(numel(events)+(1:5)),'ylim',[-3 4]);

%% Plot number of examples
figure(20); cla;
% imagesc(squeeze(mean(nExamples,2))');
% set(gca,'xtick',1:5,'ytick',1:numel(events),'yticklabel',events)
% colorbar
bar(squeeze(mean(nExamples,2)));
set(gca,'xtick',1:5)
xlabel('SqNum')
ylabel('Mean nExamples')
grid on
legend(events)
set(gcf,'Position',[13 152 576 344]);




%% CALCULATE STATS
nSubjects = numel(subjects);
nT = size(erps,2);
nChans = numel(gridOfChans);
difference = nan(nChans,nT,nTypes,nTypes,numel(events)); % sign of significant differences between classes
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
                X(:,k) = squeeze(erps(iChan,j,k,:,l));%timecourse{i,k}(:,j);
            end
            [~,~,stats] = anova1(X,types,'off');
            comparison = multcompare(stats,'alpha',0.05,'display','off','ctype','tukey-kramer');
            sig_pos = find(comparison(:,3)>0);
            sig_neg = find(comparison(:,5)<0);
            for k=1:numel(sig_pos)
                difference(i,j,comparison(sig_pos(k),1),comparison(sig_pos(k),2),l) = -1;
                difference(i,j,comparison(sig_pos(k),2),comparison(sig_pos(k),1),l) = 1;
            end
            for k=1:numel(sig_neg)
                difference(i,j,comparison(sig_neg(k),1),comparison(sig_neg(k),2),l) = 1;
                difference(i,j,comparison(sig_neg(k),2),comparison(sig_neg(k),1),l) = -1;
            end    
        end
    end    
end
close(hwait);
disp('Done!')

%% Plot Stats
tResponse = EEG_event.times;

for l=1:numel(events)
    figure(50+l); clf;
    for j=1:nTypes
        for i = 1:nChans
            subplot(nChans,nTypes,(i-1)*nTypes+j)
            cla; hold on;
            for k=1:nTypes
                plot(tResponse,difference(i,:,j,k,l)*k,'color',colors(k,:));
            end
            axis([tResponse(1), tResponse(end), -nTypes-1, nTypes+1]);
            set(gca,'xgrid','on')
            plot(get(gca,'xlim'),[0 0],'k--');
            ylabel(gridOfChans{i});
        end
        xlabel('time (ms)');
        subplot(nChans,nTypes,j);
        title(sprintf('%s vs. all others (One-Way Anova with multcompare, p=0.05)',types{j}));
        subplot(nChans,nTypes,1);
        MakeLegend(colors,types);
        MakeFigureTitle(events{l});
        set(50+l,'Position',[4 960 600 1000]);
    end
end

CascadeFigures(50+(1:numel(events)),5);

%%
for l=1:numel(events)
    figure(50+l);
    MakeFigureTitle(events{l});
    for i=1:nChans
        for j=1:nTypes
            subplot(nChans,nTypes,(i-1)*nTypes+j)
            if j==1
                ylabel(gridOfChans{i});
            else
                ylabel('');
            end
        end
    end
end