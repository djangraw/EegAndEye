% TopLevelEegGlm_script
%
% Created 4/18/13 by DJ.
% Updated 3/19/14 by DJ - tResponse in cells.


%% Set up
% % -- Squares Info
% prefix = 'sq';
% subjects = [9:11 13:15 17:19 20:27];
% glmType = 'Events-SqNum-Type-v3pt0';
% % glmType = 'LREvents-SqNum-NewType-v2pt3';
% % --Squares contrasts
% plus_events = {'aT_{0/2}', 'aT^*_{1/2}', 'aT^*_{1/2}'};
% minus_events = {'aD_{0/2}','aD_{0/2}','aT_{0/2}'};
% % plus_events = {'Compl', 'Compl', 'Compl', 'Compl', ...
% %     'New-Integ', 'Incompl'};
% % minus_events = {'New-Dist-0T','New-Dist-1T','New-Integ','Incompl',...
% %     'New-Dist-0T','New-Dist-1T'};

% % -- SquaresFix Info
% prefix = 'sf';
% subjects = [1:10 12:13];
% glmType = 'Events-SqNum-Type-v3pt1';
% % glmType = 'SqNum-Type-v2pt3';
% % --SqFix contrasts
% plus_events = {'pT_{0/2}', 'pT^*_{1/2}', 'pT^*_{1/2}'};
% minus_events = {'pD_{0/2}','pD_{0/2}','pT_{0/2}'};
% % plus_events = {'sf-Compl'};
% % minus_events = {'sf-Dist-1T'};

% -- SquaresFix3 Info
prefix = 'sf3';
subjects = 1:12;
glmType = 'Events-SqNum-Type-v3pt1';
% plus_events = {'pT_{0/3}', 'pT_{1/3}', 'pT_{1/3}', 'pT^*_{2/3}', 'pT^*_{2/3}','pT^*_{2/3}'};
% minus_events = {'pD_{0/3}','pD_{0/3}','pT_{0/3}','pD_{0/3}','pT_{0/3}','pT_{1/3}',};
plus_events = {'pD_{0/3}', 'pT_{0/3}', 'pD_{1/3}', 'pT_{1/3}', 'pD_{2/3}', 'pT^*_{2/3}'};
minus_events = repmat({''},1,numel(plus_events));

% Set up contrasts
% plus_events = {'Incompl'};
% minus_events = {'New-Dist-1T'};
% plus_events = R(1).regressor_events{end};
% minus_events = repmat({''},size(plus_events));

bandname = '';
baseline_win = []; % no baseline subtraction
multcorrect = '';

%% Load results
clear R
bandname = '';
for i=1:numel(subjects)
    fprintf('Loading subject %d/%d...\n',i,numel(subjects));
    filename = sprintf('%s-%d-GLMresults-%s%s',prefix,subjects(i),bandname,glmType);
%     filename = sprintf('%s-%d-ERPresults-%s',prefix,subjects(i),glmType);
    R(i) = load(filename);
end

%% Run Hierarchical GLM
[contrastFns, group_RF, contrastZ, group_z] = deal(cell(1,numel(plus_events)));
for i=1:numel(plus_events)
    [contrastFns{i}, group_RF{i}, contrastZ{i}, group_z{i}] = RunTopLevelGlm_Full(R,plus_events{i},minus_events{i},baseline_win,multcorrect);  
end

%% Save results
save(sprintf('GroupGlm_%s_%s_%s',prefix,glmType,datestr(now,29)),'prefix','subjects',...
    'glmType','plus_events','minus_events','baseline_win','multcorrect',...
    'contrastFns','group_RF', 'contrastZ', 'group_z');


%% v2p4


% Set up
prefix = 'sq';
subjects = [9:11 13:15 17:19];
% glmType = 'LREvents-SqNum-Type-v2pt4';
glmType = 'LREvents-Type-SqNum-v2pt4';
% prefix = 'sf';
% subjects = 1:10;
% % glmType = 'SqNum-Type-v2pt4';
% glmType = 'Type-SqNum-v2pt4';

% Set up contrasts
% plus_events = {'Incompl'};
% minus_events = {'New-Dist-1T'};
% plus_events = R(1).regressor_events{end};
% minus_events = repmat({''},size(plus_events));
% --Squares contrasts
% contrast_suffix_cell = {'CvsD1' 'CvsI','IvsD0','CvsD0','CvsInc','IncvsD1','D1vsD0'};
% contrast_events_cell = {{'Compl','New-Dist-1T'},{'Compl','New-Integ'},...
%     {'New-Integ','New-Dist-0T'},{'Compl','New-Dist-0T'},...
%     {'Compl','Incompl'},{'Incompl','New-Dist-1T'},{'New-Dist-1T','New-Dist-0T'}};
% bandnames = {'','Theta_v2','UpperAlpha_v2','Lower1Alpha_v2','Lower2Alpha_v2'};%,'Beta_v2','LowerGamma_v2'};

% --SqNum contrasts
contrast_suffix_cell = {'S1vsS5','S2vsS5','S3vsS5','S4vsS5'};
contrast_events_cell = {{'SqNum1','SqNum5'},{'SqNum2','SqNum5'},{'SqNum3','SqNum5'},{'SqNum4','SqNum5'}};
bandnames = {''};%,'Theta_v2','UpperAlpha_v2','Lower1Alpha_v2','Lower2Alpha_v2'};%,'Beta_v2','LowerGamma_v2'};



% --SqFix contrasts
% contrast_suffix_cell = {'CvsD1' 'CvsI','IvsD0','CvsD0','CvsInc','IncvsD1','D1vsD0'};
% contrast_events_cell = {{'sf-Compl','sf-Dist-1T'},{'sf-Compl','sf-Integ'},...
%     {'sf-Integ','sf-Dist-0T'},{'sf-Compl','sf-Dist-0T'},...
%     {'sf-Compl','sf-Incompl'},{'sf-Incompl','sf-Dist-1T'},{'sf-Dist-1T','sf-Dist-0T'}};
% bandnames = {'','Theta_v2','UpperAlpha_v2','Lower1Alpha_v2','Lower2Alpha_v2','Beta_v2','LowerGamma_v2'};

% --SqNum contrasts
% contrast_suffix_cell = {'S1vsS5','S2vsS5','S3vsS5','S4vsS5'};
% contrast_events_cell = {{'sf-SqNum1','sf-SqNum5'},{'sf-SqNum2','sf-SqNum5'},{'sf-SqNum3','sf-SqNum5'},{'sf-SqNum4','sf-SqNum5'}};
% bandnames = {''};%,'Theta_v2','UpperAlpha_v2','Lower1Alpha_v2','Lower2Alpha_v2'};%,'Beta_v2','LowerGamma_v2'};

% contrast_suffix_cell = {'S1vsCir','S2vsCir','S3vsCir','S4vsCir','S5vsCir'};
% contrast_events_cell = {{'sf-SqNum1','sf-Circle'},{'sf-SqNum2','sf-Circle'},{'sf-SqNum3','sf-Circle'},{'sf-SqNum4','sf-Circle'},{'sf-SqNum5','sf-Circle'}};
% bandnames = {''};%,'Theta_v2','UpperAlpha_v2','Lower1Alpha_v2','Lower2Alpha_v2'};%,'Beta_v2','LowerGamma_v2'};


baseline_win = []; % no baseline subtraction
multcorrect = '';


% Run Hierarchical GLM
% [contrastFns, group_RF, contrastZ, group_z] = deal(cell(1,numel(contrast_events_cell)));


for iBand = 1:numel(bandnames)
    fprintf('--- Band: %s ---\n',bandnames{iBand});
    bandname = bandnames{iBand};
    if ~isempty(bandname)
        bandname = strcat(bandname,'-');
    end

    if strcmp(prefix,'sf') || contrast_suffix_cell{1}(1)=='S'  
        clear R
        for j=1:numel(subjects)
            fprintf('Loading subject %d/%d...\n',j,numel(subjects));
            filename = sprintf('%s-%d-GLMresults-%s%s',prefix,subjects(j),bandname,glmType);
            foo = load(filename);
            if j==1
                R(j) = foo;
            else
                R(j) = orderfields(foo,R(1));
            end
        end
    end
    
    for i=1:numel(contrast_events_cell)

        fprintf('=== Starting %s Upper-level GLM %d/%d...\n',bandname,i,numel(contrast_events_cell));
        contrast_suffix = contrast_suffix_cell{i};
        plus_events = contrast_events_cell{i}(1); % keep in a cell
        minus_events = contrast_events_cell{i}(2); % keep in a cell

        if strcmp(prefix,'sq') && contrast_suffix_cell{1}(1)~='S'% Load results
            clear R
            for j=1:numel(subjects)
                fprintf('Loading subject %d/%d...\n',j,numel(subjects));
                filename = sprintf('%s-%d-GLMresults-%s%s%s',prefix,subjects(j),bandname,glmType,contrast_suffix);           
                foo = load(filename);
                if j==1
                    R(j) = foo;
                else
                    R(j) = orderfields(foo,R(1));
                end
            end
        end

        % Run GLM
        fprintf('=== Running %s Upper-level GLM %d/%d...\n',bandname,i,numel(contrast_events_cell));
        [contrastFns, group_RF, contrastZ, group_z] = deal(cell(1));
        [contrastFns{1}, group_RF{1}, contrastZ{1}, group_z{1}] = RunTopLevelGlm_Full(R,plus_events{1},minus_events{1},baseline_win,multcorrect);  

        % Save results
        fprintf('=== Saving %s Upper-level GLM %d/%d...\n',bandname,i,numel(contrast_events_cell));
        save(sprintf('GroupGlm_%s_%s%s%s_%s',prefix,bandname,glmType,contrast_suffix,datestr(now,29)),'prefix','subjects',...
            'glmType','plus_events','minus_events','baseline_win','multcorrect',...
            'contrastFns','group_RF', 'contrastZ', 'group_z');
    end

end

disp('Done!')







%% Make group RF scalp plots

bandnames = {''};
datestring = '2014-03-13'; % when the data was saved
switch prefix
    case 'sq'
        contrast_events_cell = {{'aT_{0/2}','aD_{0/2}'}, ...
            {'aT^*_{1/2}','aD_{0/2}'}, ...
            {'aT^*_{1/2}','aT_{0/2}'}};
        rowlabels = {'aT_{0/2} - aD_{0/2}', ...
            'aT^*_{1/2} - aD_{0/2}', ...
            'aT^*_{1/2} - aT_{0/2}'};
    case 'sf'
        contrast_events_cell = {{'pT_{0/2}','pD_{0/2}'}, ...
            {'pT^*_{1/2}','pD_{0/2}'}, ...
            {'pT^*_{1/2}','pT_{0/2}'}};
        rowlabels = {'pT_{0/2} - pD_{0/2}', ...
            'pT^*_{1/2} - pD_{0/2}', ...
            'pT^*_{1/2} - pT_{0/2}'};
    case 'sf3'
        contrast_events_cell = {{'pT_{0/3}','pD_{0/3}'}, ...
            {'pT_{1/3}','pD_{0/3}'}, ...
            {'pT_{1/3}','pT_{0/3}'},...
            {'pT^*_{2/3}','pD_{0/3}'}, ...
            {'pT^*_{2/3}','pT_{0/3}'},...
            {'pT^*_{2/3}','pT_{1/3}'}};
        rowlabels = {'pT_{0/3} - pD_{0/3}', ...
            'pT_{1/3} - pD_{0/3}', ...
            'pT_{1/3} - pT_{0/3}',...
            'pT^*_{2/3} - pD_{0/3}', ...
            'pT^*_{2/3} - pT_{0/3}',...
            'pT^*_{2/3} - pT_{1/3}'};
    otherwise
        error('prefix not recognized!');
end
contrast_suffix_cell = repmat({''},1,numel(rowlabels));    
if iscell(R(1).tResponse)
    tResponse = R(1).tResponse{R(1).iLevel};
else
    tResponse = R(1).tResponse;
end

% datestring = '2013-09-06';
% contrast_suffix_cell = {'IvsD0','CvsI','CvsD0','CvsD1','CvsInc','D1vsD0','IncvsD1'};
% rowlabels = {'ST_0_/_2 - SD_0_/_2','ST_1_/_2^* - ST_0_/_2','ST_1_/_2^* - SD_0_/_2',...
%     'ST_1_/_2^* - SD_1_/_2','ST_1_/_2^* - SD_1_/_2^*','SD_1_/_2 - SD_0_/_2',...
%     'SD_1_/_2^* - SD_1_/_2'};

% Set up
option = 'z';%'rf';%
% tPlots = 50:100:450;
tPlots = 25:50:475;
chansToPlot = {'F7','F3','FZ','F4','F8';...
                'T7','C3','CZ','C4','T8';...
                'P7','P3','PZ','P4','P8';...
                'PO7','PO3','POZ','PO4','PO8'};

for iBand = 1%1:numel(bandnames)
    fprintf('--- Band: %s ---\n',bandnames{iBand});
    bandname = bandnames{iBand};
    if ~isempty(bandname)
        bandname = strcat(bandname,'-');
    end
        
    results = nan(R(1).EEG.nbchan, numel(tPlots), numel(contrast_events_cell));
    gRF = nan(R(1).EEG.nbchan, numel(tResponse), numel(contrast_events_cell));
%     rowlabels = cell(1,numel(contrast_events_cell));
    for i=1:numel(contrast_events_cell)        
        fprintf('=== Loading %s Upper-level GLM %d/%d...\n',bandname,i,numel(contrast_events_cell));
        contrast_suffix = contrast_suffix_cell{i};
        
        load(sprintf('GroupGlm_%s_%s%s%s_%s',prefix,bandname,glmType,contrast_suffix,datestring)); 

        switch option
            case 'rf'
                %%% Response Fns
                gRF(:,:,i) = cat(3,group_RF{i});
                if isempty(bandname)
                    clim = [-2 2];
                else
                    clim = [-.2 .2];
                end
                cthresh = 0;
            case 'z'
                %%% Z-scores
                gRF(:,:,i) = cat(3,group_z{i});
                clim = [-3 3];
                cthresh = 1.96;
        end

        % Get plots
        [results(:,:,i), meanResults] = GetScalpMaps(gRF(:,:,i),tResponse,tPlots,[],cthresh);
%         rowlabels{i} = sprintf('%s - %s',plus_events{1},minus_events{1});
    end
    disp('Plotting...')
    % Make plots
%     figure%(iBand); 
    clf;
%     PlotScalpMaps(results,R(1).EEG.chanlocs,clim,tPlots,rowlabels);
    PlotScalpMaps(results,R(1).EEG.chanlocs,clim,tPlots,rowlabels);
    if strcmp(option,'rf')
        MakeFigureTitle(sprintf('Group Betas: %s, %d subj, %s%s GLM',prefix,numel(subjects),bandname,glmType));
    else
        MakeFigureTitle(sprintf('Group Z-scores: %s, %d subj, %s%s GLM',prefix,numel(subjects),bandname,glmType));
    end
    % PlotScalpMaps(cat(3,results, meanResults),R(1).EEG.chanlocs,clim,tPlots);
%     figure(numel(bandnames)+iBand); clf; 
%     PlotResponseFnsGrid(gRF,rowlabels,tResponse,chanlocs,chansToPlot);
%     if strcmp(option,'rf')
%         MakeFigureTitle(sprintf('Group Betas: %s, %d subj, %s%s GLM',prefix,numel(subjects),bandname,glmType));
%     else
%         MakeFigureTitle(sprintf('Group Z-scores: %s, %d subj, %s%s GLM',prefix,numel(subjects),bandname,glmType));
%     end
    
    disp('Done!')
%     pause

end

%% Annotate figure
M = length(tPlots);
N = size(gRF,3);
for i=1:M
    for j=1:N
        subplot(N+1,M,(j-1)*M+i);
        title(sprintf('%s - %s, %.1f ms',plus_events{j},minus_events{j},tPlots(i)));
    end
    subplot(N+1,M,N*M+i);
    title(sprintf('Mean, %.1f ms',tPlots(i)));
end
MakeFigureTitle(sprintf('Group Betas: %s, %d subj, %s GLM',prefix,numel(subjects),glmType));

%% Make z score movie
% Make movie
i=1;
MakeTopoMovie(group_z{i},tResponse,R(1).EEG.chanlocs,[],[-3 3]);
% Annotate
MakeFigureTitle(sprintf('%s GLM, %s - %s contrast (mc correction: %s)',glmType,plus_events,minus_events,multcorrect),0);
