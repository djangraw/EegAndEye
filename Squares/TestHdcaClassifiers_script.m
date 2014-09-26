% TestHdcaClassifiers_script
%
% Created 3/17/14 by DJ.

% R_sf3_type = LoadAllGlmResults('sf3','Events-SqNum-Type-v3pt1');

% Set up
train_prefix = 'sf3';%'sf'; %
train_events = {'pT_{0/3}','pT^*_{2/3}'};% {'pT_{0/2}','pT^*_{1/2}'}; %
test_events = {'pD_{2/3}','pD_{1/3}','pD_{0/3}','pT_{0/3}','pT_{1/3}','pT^*_{2/3}'};% {'pT_{0/2}','pT^*_{1/2}'};%
eval(sprintf('R_train = R_%s_type;',train_prefix));
% twl_ms = 50;
% two_ms = 0:20:700;
twl_ms = 150;
two_ms = 200;

%% Run Classifier
truth = deal(cell(1,numel(R_train)));
[y,w,v,fwdModel,y_level1] = deal(cell(1,numel(R_train)));

for i=1:numel(R_train)
%     if ~isempty(y{i}), continue; end
    fprintf('--- Subject %d/%d...\n',i,numel(R_train));
    % Subtract out other factors
    EEG = R_train(i).EEG;
    EEG = SubtractOutGlmResponses(EEG,R_train(i).responseFns{1},R_train(i).influence{1},R_train(i).regressor_events{1});
    EEG = SubtractOutGlmResponses(EEG,R_train(i).responseFns{2},R_train(i).influence{1},R_train(i).regressor_events{2});  
    
    % Extract epochs   
    EEG0 = pop_epoch( EEG, train_events(1), [-0.5 1], 'newname', train_events{1}, 'epochinfo', 'yes');    
    EEG1 = pop_epoch( EEG, train_events(2), [-0.5 1], 'newname', train_events{2}, 'epochinfo', 'yes');
    
    % Convert twl and two to indices
    t_ms = EEG0.times;
    twl = twl_ms/1000*EEG0.srate; % get window width in samples
    two = round(interp1(t_ms,1:length(t_ms),two_ms)); % get nearest indices
    
    % Assemble data & truth vectors
    data = cat(3,EEG0.data,EEG1.data);
    truth{i} = [zeros(EEG0.trials,1); ones(EEG1.trials,1)]';
    
    % Run classifier
    [y{i},w{i},v{i},fwdModel{i},y_level1{i}] = RunHybridHdcaClassifier_LR(data,truth{i},twl,two,'10fold');
%     [y2{i},w2{i},v2{i},fwdModel2{i},y_level12{i}] = RunHybridHdcaClassifier(data,truth{i},twl,two,'10fold');
    
    figure;
    PlotHybridHdcaClassifier(fwdModel{i},v{i},EEG0.chanlocs,two_ms + twl_ms/2);
end
disp('--- Done!')


%% Plot averages

allFwdModels = squeeze(mean(cat(4,fwdModel{:}),3)); % allFwdModels(:,:,i) is subject i's mean fm across folds.
allVs = permute(mean(cat(4,v{:}),3),[1 2 4 3]); % allVs(:,:,i) is subject i's mean v across folds.
%
figure;
PlotHybridHdcaClassifier(allFwdModels,allVs,EEG0.chanlocs,two_ms + twl_ms/2);

%% Plot Az values
Az_level1 = nan(numel(R_train),numel(two_ms));
Az_level2 = nan(numel(R_train),1);
for i=1:numel(R_train)
    for j=1:numel(two_ms)
        Az_level1(i,j) = rocarea(y_level1{i}(:,j),truth{i});
    end
    Az_level2(i) = rocarea(y{i},truth{i});
end

% Plot Az's
figure;
errorbar(two_ms + twl_ms/2, mean(Az_level1,1),std(Az_level1,[],1)/sqrt(numel(R_train)));
% Annotate
set(gca,'xgrid','on')
ylim([.45 .7])
hold on
plot(get(gca,'xlim'),[.5 .5],'k--');
% Label
xlabel('Time of bin center (ms)');
ylabel('Single-bin AUC');
title(sprintf('Multi-bin AUC is %.3f +/- %.3f\n(mean +/- stderr across %d subjects)',mean(Az_level2),std(Az_level2)/sqrt(numel(R_train)),numel(R_train)));

%% Show time course of T1, T2, T3

figure;
nRows = ceil(numel(R_train)/2);
nCols = 4;
iALLOWED_BINS = 1;%10:15; % 4:7; % 1:8;%
[wMax,fmMax] = deal(cell(1,numel(R_train)));
nEvents = numel(test_events);
y_cell = cell(nEvents,numel(R_train));
method = 'rf';

for i=1:numel(R_train)
    fprintf('--- Subject %d/%d...\n',i,numel(R_train));
    
    % Get weights and fwd model from best window
    [~,iABMax] = max(Az_level1(i,iALLOWED_BINS)); % best weights
    iMax = iALLOWED_BINS(iABMax);
    wMax{i} = mean(mean(w{i}(:,iMax),3),4); % mean across folds for best weights    
    fmMax{i} = mean(fwdModel{i}(:,iMax),3); % mean across folds for best weights
    wMax{i} = wMax{i}/sqrt(sum(wMax{i}.^2)); % NORMALIZE!
    
    switch method
        case 'erp'
            % Subtract out other factors
            EEG = R_train(i).EEG;
            EEG = SubtractOutGlmResponses(EEG,R_train(i).responseFns{1},R_train(i).influence{1},R_train(i).regressor_events{1});
            EEG = SubtractOutGlmResponses(EEG,R_train(i).responseFns{2},R_train(i).influence{2},R_train(i).regressor_events{2});  

            for j=1:nEvents                
                % Extract epochs
                EEG0 = pop_epoch( EEG, test_events(j), [-0.5 1], 'newname', test_events{j}, 'epochinfo', 'yes');
                % Apply weights of interest to data
                y_cell{j,i} = wMax{i}'*mean(EEG0.data,3);
            end
            t_ms = EEG0.times;    
                

        case 'rf'
                        
            for j=1:nEvents
                i0 = find(strcmp(test_events{j},R_train(i).regressor_events{3}));
                y_cell{j,i} = wMax{i}'*R_train(i).responseFns{3}(:,:,i0);
            end
            t_ms = R_train(i).tResponse{3};
            
    end
    % Plot scalp map
    subplot(nRows,nCols,i*2-1);
    topoplot(fmMax{i},EEG0.chanlocs);
    title(sprintf('Subject %d: %d ms',i,two_ms(iMax) + twl_ms/2));
    
    % Plot timecourse
    subplot(nRows,nCols,i*2);
    cla;    
    plot(t_ms,cat(1,y_cell{:,i}));
    hold on;
    plot([t_ms(1) t_ms(end)],[0 0],'k--');
    PlotVerticalLines(two_ms(iMax) + [0 twl_ms],'m:');
    title(sprintf('Max Az = %.2f',Az_level1(i,iMax)));
    xlabel('time (ms)');
    ylabel('amplitude (uV)');
    xlim([t_ms(1) t_ms(end)]);

end

legend(test_events);

% Plot mean across subjects
figure;

clf;
subplot(1,2,1);
topoplot(mean(cat(3,fmMax{:}),3),EEG0.chanlocs);
title(sprintf('Mean forward model across %d %s subjects',numel(R_train),train_prefix));
subplot(1,2,2);
ymat = [];
for j=1:nEvents
    ymat = cat(3,ymat,cat(4,y_cell{j,:})); % dim3 = events, dim4 = subj
end
Cmap = GetSquaresEventColormap(test_events);
PlotResponseFns(ymat,test_events,t_ms,struct('labels',''),'',Cmap)

xlim([0 750]);
title('Mean timecourse with individual-subject weights');
xlabel('time (ms)');
ylabel('amplitude (uV)');

% Perform stats check on y2>y1>y0
subplot(1,2,2); hold on; % superimpose significance stars on plot from last code section
ymean = mean(ymat,4);
[h,p] = deal(nan(nEvents-1,length(t_ms)));
iStars = [1:nEvents-1];
for j=1:nEvents-1
    for i=1:length(t_ms)
%             [h(1,i),p(1,i)] = ttest(squeeze(ymat(1,i,2,:)), squeeze(ymat(1,i,1,:)),.05,'right');
        [p(j,i),h(j,i)] = signrank(squeeze(ymat(1,i,j+1,:)), squeeze(ymat(1,i,j,:)),'tail','right');
    end
    % Plot stars  
    plot(t_ms(h(j,:)==1),ymean(1,h(j,:)==1,iStars(j)),'*','Color',Cmap(iStars(j),:),'linewidth',3);
end




%% Apply mean weights to other datasets
% fmMean = mean(cat(3,fmMax{:}),3);
wMean = mean(cat(3,wMax{:}),3);
prefixes = {'sf3','sf','sq'};
target_events = {{'pD_{2/3}','pD_{1/3}','pD_{0/3}','pT_{0/3}','pT_{1/3}','pT^*_{2/3}'},...
    {'pD_{1/2}','pD_{0/2}','pT_{0/2}','pT^*_{1/2}'},...
    {'aD_{1/2}','aD_{0/2}','aT_{0/2}','aT^*_{1/2}'}};
smoothing_sigma = 5;
smoothing_type = 'full'; % 'half';
% target_events = {{'pT_{0/3}','pT_{1/3}','pT^*_{2/3}'},...
%     {'pT_{0/2}','pT^*_{1/2}'},...
%     {'aT_{0/2}','aT^*_{1/2}'}};
method = 'rf';

for iExp=1:numel(prefixes)
    eval(sprintf('R = R_%s_type;',prefixes{iExp}));
    figure;
    nRows = ceil(sqrt(numel(R)));
    nCols = ceil(numel(R)/nRows);    
    legendstr = target_events{iExp};
    nEvents = numel(legendstr);
    y_cell = cell(nEvents,numel(R));
    
    for i=1:numel(R)
        fprintf('--- Subject %d/%d...\n',i,numel(R));
        switch method
            case 'erp'
                % Subtract out other factors
                EEG = R(i).EEG;
                EEG = SubtractOutGlmResponses(EEG,R(i).responseFns{1},R(i).influence{1},R(i).regressor_events{1});
                EEG = SubtractOutGlmResponses(EEG,R(i).responseFns{2},R(i).influence{2},R(i).regressor_events{2});  

                % Extract epochs
                for j=1:nEvents
                    EEG0 = pop_epoch( EEG, legendstr(j), [-0.5 1], 'newname', legendstr{j}, 'epochinfo', 'yes');
                    y_cell{j,i} = wMean'*mean(EEG0.data,3);
                end
            case 'rf'
                for j=1:nEvents
                    % Find events of interest
                    i0 = find(strcmp(legendstr{j},R(i).regressor_events{3}));
                    % Apply weights of interest to data
                    y_cell{j,i} = wMean'*R(i).responseFns{3}(:,:,i0);
                
                end                
        end

        % Plot timecourse
        subplot(nRows,nCols,i);
        cla;    
        plot(t_ms,cat(1,y_cell{:,i}));
        hold on;
        plot([t_ms(1) t_ms(end)],[0 0],'k--');
        xlabel('time (ms)');
        ylabel('amplitude (uV)');
    
    end
    legend(legendstr);
    
    
    % Plot mean across subjects
    figure; clf;
    subplot(1,2,1);
    topoplot(mean(cat(3,fmMax{:}),3),EEG0.chanlocs);
    title(sprintf('Mean forward model across %d %s subjects',numel(R_train),train_prefix));
    subplot(1,2,2);
    ymat = [];
    for j=1:nEvents
        ymat = cat(3,ymat,cat(4,y_cell{j,:})); % dim3 = events, dim4 = subj
    end
    
    Cmap = GetSquaresEventColormap(legendstr);
%     PlotResponseFns(ymat,legendstr,t_ms,struct('labels',''),'',Cmap)
    ymat_smooth = SmoothData(ymat,smoothing_sigma,smoothing_type);
    PlotResponseFns(ymat_smooth,legendstr,t_ms,struct('labels',''),'',Cmap)
    % plot(t_ms,meanYmat);
    xlim([0 750]);
    title(sprintf('Mean timecourse across %d %s subjects',numel(R),prefixes{iExp}));
    xlabel('time (ms)');
    ylabel('amplitude (uV)');
    MakeLegend(Cmap,legendstr,2);
    
    
    
    % Perform stats check on y2>y1>y0
%     ymean = mean(ymat,4);
    ymean_smooth = mean(ymat_smooth,4);
    [h,p] = deal(nan(nEvents-1,length(t_ms)));
    subplot(1,2,2); hold on; % superimpose significance stars on plot from last code section
    iStars = [1:nEvents-1]; %[1 3];
    for j=1:nEvents-1
        for i=1:length(t_ms)
%             [h(1,i),p(1,i)] = ttest(squeeze(ymat(1,i,2,:)), squeeze(ymat(1,i,1,:)),.05,'right');
            [p(j,i),h(j,i)] = signrank(squeeze(ymat(1,i,j+1,:)), squeeze(ymat(1,i,j,:)),'tail','right');
        end
        % Plot stars  
%         plot(t_ms(h(j,:)==1),ymean(1,h(j,:)==1,iStars(j)),'*','Color',Cmap(iStars(j),:),'linewidth',3);
        plot(t_ms(h(j,:)==1),ymean_smooth(1,h(j,:)==1,iStars(j)),'*','Color',Cmap(iStars(j),:),'linewidth',3);
    end
    

    
    
end

%% Apply mean weights to other datasets - ONE FIGURE
% fmMean = mean(cat(3,fmMax{:}),3);
wMean = mean(cat(3,wMax{:}),3);
prefixes = {'sf3','sf','sq'};
% target_events = {{'pD_{0/3}','pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pT^*_{2/3}'},...
%     {'pD_{0/2}','pT_{0/2}','pD_{1/2}','pT^*_{1/2}'},...
%     {'aD_{0/2}','aT_{0/2}','aD_{1/2}','aT^*_{1/2}'}};

smoothing_sigma = 2.5;
smoothing_type = 'full'; % 'half';
% target_events = {{'pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pT_{+/3}'},...
%     {'pT_{0/2}','pT^*_{1/2}','pT_{+/2}'},...
%     {'aT_{0/2}','aT^*_{1/2}','aT_{+/2}'}};
target_events = {{'pD_{2/3}','pD_{+/3}','pT^*_{2/3}','pT_{+/3}'},...
    {'pD_{1/2}','pD_{+/2}','pT^*_{1/2}','pT_{+/2}'},...
    {'aD_{1/2}','aD_{+/2}','aT^*_{1/2}','aT_{+/2}'}};

method = 'rf';

figure;
set(gcf,'Position',[1212        1283         380         204]);
topoplot(mean(cat(3,fmMax{:}),3),EEG0.chanlocs);

figure;

for iExp=1:numel(prefixes)
    eval(sprintf('R = R_%s_type;',prefixes{iExp}));
%     figure;
    nRows = ceil(sqrt(numel(R)));
    nCols = ceil(numel(R)/nRows);    
    legendstr = target_events{iExp};
    nEvents = numel(legendstr);
    y_cell = cell(nEvents,numel(R));
    
    for i=1:numel(R)
        fprintf('--- Subject %d/%d...\n',i,numel(R));
        switch method
            case 'erp'
                % Subtract out other factors
                EEG = R(i).EEG;
                EEG = SubtractOutGlmResponses(EEG,R(i).responseFns{1},R(i).influence{1},R(i).regressor_events{1});
                EEG = SubtractOutGlmResponses(EEG,R(i).responseFns{2},R(i).influence{2},R(i).regressor_events{2});  

                % Extract epochs
                for j=1:nEvents
                    EEG0 = pop_epoch( EEG, legendstr(j), [-0.5 1], 'newname', legendstr{j}, 'epochinfo', 'yes');
                    y_cell{j,i} = wMean'*mean(EEG0.data,3);
                end
            case 'rf'
                for j=1:nEvents
                    % Find events of interest
                    i0 = find(strcmp(legendstr{j},R(i).regressor_events{3}));
                    % Apply weights of interest to data
                    y_cell{j,i} = wMean'*R(i).responseFns{3}(:,:,i0);
                
                end                
        end

    
    end

    ymat = [];
    for j=1:nEvents
        ymat = cat(3,ymat,cat(4,y_cell{j,:})); % dim3 = events, dim4 = subj
    end
    
    Cmap = GetSquaresEventColormap(legendstr);
%     PlotResponseFns(ymat,legendstr,t_ms,struct('labels',''),'',Cmap)
    ymat_smooth = SmoothData(ymat,smoothing_sigma,smoothing_type);
    subplot(numel(prefixes),1,iExp);
    PlotResponseFns(ymat_smooth,legendstr,t_ms,struct('labels',''),'',Cmap)
    % plot(t_ms,meanYmat);
    xlim([0 750]);
    title(sprintf('Mean timecourse across %d %s subjects',numel(R),prefixes{iExp}));
    xlabel('time (ms)');
    ylabel('amplitude (uV)');
%     MakeLegend(Cmap,legendstr,2);
    
    
    
    % Perform stats check on y2>y1>y0
%     ymean = mean(ymat,4);
%     hold on;
%     ymean_smooth = mean(ymat_smooth,4);
%     [h,p] = deal(nan(nEvents-1,length(t_ms)));
% %     subplot(1,2,2); hold on; % superimpose significance stars on plot from last code section
%     iStars = [1:nEvents-1]; %[1 3];
%     for j=1:nEvents-1
%         for i=1:length(t_ms)
% %             [h(1,i),p(1,i)] = ttest(squeeze(ymat(1,i,2,:)), squeeze(ymat(1,i,1,:)),.05,'right');
%             [p(j,i),h(j,i)] = signrank(squeeze(ymat(1,i,j+1,:)), squeeze(ymat(1,i,j,:)),'tail','right');
%         end
%         % Plot stars  
% %         plot(t_ms(h(j,:)==1),ymean(1,h(j,:)==1,iStars(j)),'*','Color',Cmap(iStars(j),:),'linewidth',3);
%         plot(t_ms(h(j,:)==1),ymean_smooth(1,h(j,:)==1,iStars(j)),'*','Color',Cmap(iStars(j),:),'linewidth',3);
%     end
    

    
    
end

% set(gcf,'Position',[1554         536         872         960]);
set(gcf,'Position',[1899         917         578         586]);
subplot(3,1,1);
title('');
xlabel('')
set(gca,'xticklabel',{})
ylabel('Passive-3');
% ylim([-3 3])

subplot(3,1,2);
title('');
xlabel('')
set(gca,'xticklabel',{})
ylabel('Passive-2');

subplot(3,1,3);
title('');
ylabel('Active-2');


%% Stats

iTimes = find(t_ms>275 & t_ms<325);
[p,h] = deal(cell(1,3));
for iExp = 1:3
    eval(sprintf('R = R_%s_type;',prefixes{iExp}));
    legendstr = target_events{iExp};
    nEvents = numel(legendstr);
    Cmap = GetSquaresEventColormap(legendstr);
    y_cell = cell(nEvents,numel(R));
    
    for i=1:numel(R)
        fprintf('--- Subject %d/%d...\n',i,numel(R));

        for j=1:nEvents
            % Find events of interest
            i0 = find(strcmp(legendstr{j},R(i).regressor_events{3}));
            % Apply weights of interest to data
            y_cell{j,i} = wMean'*R(i).responseFns{3}(:,:,i0);

        end                

    
    end

    
    ymat = [];
    for j=1:nEvents
        ymat = cat(3,ymat,cat(4,y_cell{j,:})); % dim3 = events, dim4 = subj
    end

    [p{iExp},h{iExp}] = deal(nan(1,size(ymat,3)-1));
    for j=1:size(ymat,3)-1
        [p{iExp}(j),h{iExp}(j)] = signrank(squeeze(mean(ymat(1,iTimes,j+1,:),2)), squeeze(mean(ymat(1,iTimes,j,:),2)));%,'tail','right');
    end
    
    subplot(3,1,iExp); cla;
    hold on;
    for j=1:size(ymat,3)
        bar(j,mean(mean(ymat(1,iTimes,j,:))),'facecolor',Cmap(j,:));
    end
    nSubj = size(ymat,4);
    errorbar(mean(mean(ymat(1,iTimes,:,:),2),4),std(mean(ymat(1,iTimes,:,:),2),[],4)/sqrt(nSubj),'k.');
    set(gca,'xtick',1:size(ymat,3),'xticklabel',legendstr)
    xlabel('Event Types')
    ylabel('Amplitude (uV)');
%     ylim([400 520])

    % Display anovan result
    X = squeeze(mean(ymat(1,iTimes,:,:),2));
    grp = meshgrid(1:size(X,1),1:size(X,2))';
    Y = reshape(X',numel(X),1);
    GROUP = reshape(grp',numel(grp),1);
    [pN{iExp},tN{iExp},statsN{iExp}] = anovan(Y,GROUP);
    fprintf('%s: p = %g\n',prefixes{iExp},pN{iExp})
end
% Resize figure
set(gcf,'Position',[1899         917         400         586]);
% Annotate figure
subplot(3,1,1)
ylim([-10 10]);
limits = get(gca,'xlim');
subplot(3,1,2);
ylim([-5 5]);
xlim(limits);
subplot(3,1,3);
xlim(limits);
ylim([-2 2]);

%% Get timing matches betwen extra targets and template
for iExp=1:numel(prefixes)
    eval(sprintf('R = R_%s_type;',prefixes{iExp}));
%     figure;
    nRows = ceil(sqrt(numel(R)));
    nCols = ceil(numel(R)/nRows);    
    legendstr = target_events{iExp};
    nEvents = numel(legendstr);
%     x_cell = cell(nEvents,numel(R));
    x_mat = zeros([size(R(1).responseFns{3},1),size(R(1).responseFns{3},2), nEvents,numel(R)]);
%     x_mat = zeros([1,size(R(1).responseFns{3},2), nEvents,numel(R)]);
    for i=1:numel(R)
        fprintf('--- Subject %d/%d...\n',i,numel(R));
        for j=1:nEvents
            % Find events of interest
            i0 = find(strcmp(legendstr{j},R(i).regressor_events{3}));
            % Apply weights of interest to data
%             x_cell{j,i} = R(i).responseFns{3}(:,:,i0);
            x_mat(:,:,j,i) = R(i).responseFns{3}(:,:,i0);
%             x_mat(:,:,j,i) = wMean'*R(i).responseFns{3}(:,:,i0);
        end   
    end
    
    matchStrength = [];
    iTimes = 31:51;
    for i=1:size(x_mat,4) % subjects
%         matchStrength(i,:,1) = UpdateTemplateMatchStrength(x_mat(:,:,1,i),x_mat(:,iTimes,1,i));
        matchStrength(i,:,1) = UpdateTemplateMatchStrength(x_mat(:,:,2,i),x_mat(:,iTimes,1,i));
        matchStrength(i,:,2) = UpdateTemplateMatchStrength(x_mat(:,:,4,i),x_mat(:,iTimes,3,i));
    end
    figure(1);
    subplot(3,1,iExp); cla;
    nTimes = length(iTimes);
    imagesc(t_ms(1:(end-nTimes))-t_ms(iTimes(1)),1:numel(R),matchStrength(:,:,2));
    
    figure(2); 
    subplot(3,1,iExp); cla;
    plot(t_ms(1:(end-nTimes))-t_ms(iTimes(1)),squeeze(mean(matchStrength,1)),'linewidth',2);   
    hold on
    plot([0 0],get(gca,'ylim'),'k','linewidth',2)
    ylabel(prefixes{iExp})
    xlabel('offset (ms)')
    xlim([-240 240])
%     legend('Distractors','Targets')
end
subplot(3,1,1);
xlabel('');
set(gca,'xticklabel',{})
subplot(3,1,2);
xlabel('');
set(gca,'xticklabel',{})

set(gcf,'Position',[1899         917         400         586]);
