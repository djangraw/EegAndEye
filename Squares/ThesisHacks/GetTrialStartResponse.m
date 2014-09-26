% function GetTrialStartResponse
% Get TrialStart Response
% Created 4/19/14 by DJ.

%% Apply mean weights to other datasets - ONE FIGURE

prefixes = {'sq','sf','sf3'};
nuisance_events = {{'Square','Circle'},...
    {'sf-Square','sf-Circle'},...
    {'sf3-Square','sf3-Circle'}};
target_events = {{'TrialStart-T','TrialStart-D'},...
    {'TrialStart-T','TrialStart-D'},...
    {'TrialStart-T','TrialStart-D'}};

%%
y_mat = cell(1,numel(prefixes));
for iExp=1:numel(prefixes)
    eval(sprintf('R = R_%s_type;',prefixes{iExp}));
    R = UpdateGlmResultsFormat(R);
%     figure;
    nRows = ceil(sqrt(numel(R)));
    nCols = ceil(numel(R)/nRows);    
    legendstr = target_events{iExp};
    nEvents = numel(legendstr);
    y_cell = cell(1,numel(R));
    
    for i=1:numel(R)
        fprintf('--- %s Subject %d/%d...\n',prefixes{iExp},i,numel(R));
        % Subtract out other factors
        EEG = R(i).EEG;
        EEG = SubtractOutGlmResponses(EEG,R(i).responseFns{1}(:,:,1:2),R(i).influence{1},nuisance_events{iExp});
        
        % old way: square onset for all three
%         EEG0 = pop_epoch( EEG, legendstr, [-1 4], 'newname', 'foo', 'epochinfo', 'yes');
        % new way: when subject knows trial has started
        switch prefixes{iExp}
            case 'sq'
                EEG0 = pop_epoch( EEG, legendstr, [-1 4], 'newname', 'foo', 'epochinfo', 'yes');
                t_ms = EEG0.times;
            otherwise
                EEG0 = pop_epoch( EEG, legendstr, [-1 4]-0.2, 'newname', 'foo', 'epochinfo', 'yes'); % offset of fix cross
        end
        y_cell{i} = mean(EEG0.data,3);
        
    end

    y_mat{iExp} = cat(3,y_cell{:}); % dim3 = subj
   
end

%% Combine
group_RF_all = nan(size(y_mat{1},1),size(y_mat{1},2),3);
group_SEM_all = nan(size(y_mat{1},1),size(y_mat{1},2),3);
for iExp=1:3
    group_RF_all(:,:,iExp) = mean(y_mat{iExp},3);
    group_SEM_all(:,:,iExp) = std(y_mat{iExp},[],3)/sqrt(size(y_mat{iExp},3));
end

chanlocs = EEG.chanlocs;
tResponse = t_ms;

%% Plot ERPs
chansToPlot = {'FZ';'CZ';'PZ';'OZ'};
colors = {'r','g','b'};
legendstr = {'Active-2','Passive-2','Passive-3'};

% figure;
clf;
set(gcf,'Position',[0 623 704 882]);

PlotResponseFnsGrid(group_RF_all,legendstr,tResponse,chanlocs,chansToPlot,colors);
for i=1:4
    subplot(4,1,i);
    ylabel(chansToPlot{i});
    xlim([-500 1000]);
    ylim([-2 2]) % to match target plot
    title('');
    set(gca,'xtick',-500:100:1000);
    if i<4
        xlabel('');
        set(gca,'xticklabel',{});
    end
end
%% Plot Scalp Maps
tBinCenters = [25 250 375 550 875];%-450:100:1000;
% tBinCenters = [-25 175 250 350 675];%-450:100:1000;%[-25 175 275 675];%37.5:75:500;
tBinWidth = 50;%50;%75;
clim = [-2.65 2.65]; % to match target plot

% figure;
clf;
% set(gcf,'Position',[684 1001 447 428]); % for 3 times
set(gcf,'Position',[684 1001 780 428]); % for 5 times

[sm_all,sm] = GetScalpMaps(group_RF_all,tResponse,tBinCenters,tBinWidth);
PlotScalpMaps(sm_all,chanlocs,clim,tBinCenters-tBinWidth/2,legendstr);

%% Find times of next events
next_events = {'TrialStart','sf-SqNum2','sf-SqNum2'};
for iExp = 1:3
    eval(sprintf('R = R_%s_type;',prefixes{iExp}));
    R = UpdateGlmResultsFormat(R);
    for i=1:numel(R)
        fprintf('--- %s Subject %d/%d...\n',prefixes{iExp},i,numel(R));
        EEG = R(i).EEG;
        tNext{iExp}{i} = nan(1,EEG.trials);
        if iExp==1
            tSaccStart{i} = nan(1,EEG.trials);
        end
        for j=1:EEG.trials
            iNext = find(strcmp(next_events{iExp},EEG.epoch(j).eventtype),1);
            if ~isempty(iNext)
                tNext{iExp}{i}(j) = EEG.epoch(j).eventlatency{iNext};
                if iExp==1
                    iSaccStart = find(strcmp('FSL',EEG.epoch(j).eventtype(1:iNext)) | strcmp('FSR',EEG.epoch(j).eventtype(1:iNext)),1,'first');
                    if ~isempty(iSaccStart)
                        tSaccStart{i}(j) = EEG.epoch(j).eventlatency{iSaccStart};
                    end
                end
            end
        end
    end
end

%% histogram them
tHist = 0:10:1000;
yHist = zeros(3,length(tHist));
for iExp = 1:3
    tMedian = zeros(1,length(tNext{iExp}));
    for i=1:length(tNext{iExp})
        yHist(iExp,:) = hist(tNext{iExp}{i},tHist);
        tMedian(i) = nanmedian(tNext{iExp}{i});
    end
    tMedian_all(iExp) = mean(tMedian);
    tSEM_all(iExp) = std(tMedian)/sqrt(length(tMedian));
end

    tSSMedian = zeros(1,length(tSaccStart));
    for i=1:length(tSaccStart)
        tSSMedian(i) = nanmedian(tSaccStart{i});
    end
    tSSMedian_all = mean(tSSMedian);
    tSSSEM_all = std(tSSMedian)/sqrt(length(tSSMedian));


figure;
subplot(2,1,1);
plot(tHist,yHist);
subplot(2,1,2);
errorbar(tMedian_all,tSEM_all);
fprintf('Median across trials, mean across subjects: %g +/- %g\n',tMedian_all(1),tSEM_all(1));
fprintf('(SaccStart median across trials, mean across subjects: %g +/- %g\n',tSSMedian_all(1),tSSSEM_all(1));