% TEMP_GetAllSvdResults_script
%
% Created 10/2/13 by DJ for one-time use.

%% Set up
nComponents = 3;
nEventTypes = 6;
timeWindow = [0 500];
componentSwap = [];
doPlot = 1;

%% Get


% prefix = 'sq';
% subjects = [9:11 13:15 17:19 20:27];
% % glmTypes = {'LREvents-SqNum-Type-v2pt4'};
% % glmTypes = {'LREvents-Type-SqNum-v2pt4'};
% glmTypes = {'Events-SqNum-Type-v3pt0'};

% prefix = 'sf';
% subjects = [1:10 12:13];
% % glmTypes = {'SqNum-Type-v2pt4'};
% glmTypes = {'Events-SqNum-Type-v3pt1'};

prefix = 'sf3';
subjects = 1:12;
glmTypes = {'Events-SqNum-Type-v3pt1'};

% bandnames = {'', 'theta_v2' 'lower1alpha_v2' 'lower2alpha_v2' 'upperalpha_v2','delta_v2' 'beta_v2' 'lowergamma_v2'};
% bandnames = {'', 'Theta_v2-' 'Lower1Alpha_v2-' 'Lower2Alpha_v2-' 'UpperAlpha_v2-','Delta_v2-' 'Beta_v2-' 'LowerGamma_v2-'};
bandnames = {''};

[weights, timecourse, avgW, avgTC, chanlocs,tResponse,steW,steTC,eventnames] = deal(cell(numel(glmTypes),numel(bandnames)));
for iGlmType = 1:numel(glmTypes)
    glmType = glmTypes{iGlmType};
    for iBand = 1:numel(bandnames);
        
        bandname = bandnames{iBand};
        fprintf('-----Running SVD for %s glmType %s, band %s\n',prefix, glmType,bandname);
        
%         if ~exist(sprintf('%s-%d-GLMresults-%s%s.mat',prefix,subjects(1),bandname,glmType),'file')
%             disp('SKIPPING!');
%             continue;
%         end
        
        % Load
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
        % Run SVD
        disp('Running SVD...')
        [weights{iGlmType,iBand}, timecourse{iGlmType,iBand}, ...
            avgW{iGlmType,iBand}, avgTC{iGlmType,iBand}, chanlocs{iGlmType,iBand}, ...
            tResponse{iGlmType,iBand}, steW{iGlmType,iBand}, ...
            steTC{iGlmType,iBand},eventnames{iGlmType,iBand}] = ...
            GetGroupSvdResults(R,nComponents,nEventTypes,timeWindow,componentSwap,doPlot);
%         % Plot SVD results
%         if ~doPlot
%             figure;
%             PlotGroupSvdResults(avgW{iGlmType,iBand}, avgTC{iGlmType,iBand}, ...
%                 chanlocs{iGlmType,iBand}, tResponse{iGlmType,iBand}, ...
%                 steW{iGlmType,iBand}, steTC{iGlmType,iBand}, eventnames{iGlmType,iBand});
%             MakeFigureTitle(sprintf('Figure %d: %s-%s%s',gcf,prefix,bandname,glmType),0);
%         end
        
        disp('Done!');
    end
end
%% Save
% filename = 'TEMP_SVDresults_sf-sqnum-v2pt4';
filename = sprintf('TEMP_SVDresults_%s-type-v3pt1',prefix);
save(filename,'prefix','subjects','glmTypes','bandnames','weights', 'timecourse', 'avgW', 'avgTC', 'chanlocs','tResponse','steW','steTC','eventnames');

%% Plot

for iGlmType = 1:numel(glmTypes)
    glmType = glmTypes{iGlmType};
% iFig = 7;
    for iBand = 1%:numel(bandnames);

        bandname = bandnames{iBand};
        fprintf('-----Plotting SVD for %s glmType %s, band %s\n',prefix, glmType,bandname);

%         if ~exist(sprintf('%s-%d-GLMresults-%s%s.mat',prefix,subjects(1),bandname,glmType),'file')
%             disp('SKIPPING!');
%             continue;
%         end

        if strcmp(glmType(end-9:end-6),'Type')  
            legendnames = R(1).regressor_events{R(1).iLevel}(1:nEventTypes);
            Cmap = zeros(length(legendnames),3);
            for j = 1:length(legendnames)
                switch legendnames{j}(2:end-2) % exclude leading a/p and trailing 2/3
                    case 'D_{0/'
                        Cmap(j,:) = [0 0 1]; %B
%                         Cmap(j,:) = [0 .75 0]; % G
                    case 'D_{1/'
                        Cmap(j,:) = [0 .75 .5]; %BG
%                         Cmap(j,:) = [0 0 1]; % B
                    case {'D_{2/' 'D_{-/'}
                        Cmap(j,:) = [0 1 1]; %C
%                         Cmap(j,:) = [.5 0 .75]; % V
                    case 'D^*_{-/' 
                        Cmap(j,:) = [0 .75 0]; %G
%                         Cmap(j,:) = [0 1 1]; % C
                    case 'D_{+/'                        
                        Cmap(j,:) = [1 1 1]*.75; % LG                    
                    case 'T_{0/'
                        Cmap(j,:) = [1 0 0]; %R
%                         Cmap(j,:) = [1 1 0]; % Y
                    case {'T_{1/' 'T^*_{1/'}
                        Cmap(j,:) = [.5 0 .75]; %V
%                         Cmap(j,:) = [1 .5 0]; % O
                    case {'T^*_{2/'}
                        Cmap(j,:) = [1 0 1]; % M
%                         Cmap(j,:) = [1 0 0]; % R
                    case 'T_{+/'
                        Cmap(j,:) = [1 1 0]; % K
%                         Cmap(j,:) = [1 0 1]; % M
                    case 'T_{-/'
                        Cmap(j,:) = [1 1 1]*.25; % DG                    
                    otherwise
                        Cmap(j,:) = [j j j]/length(legendnames); % gray
                end
            end
                        
%             if strcmp(prefix,'sf')
%                 Cmap = [0 0 1; 0 1 1; 1 0 0; 1 0 1; 0 0.75 0];
%                 legendnames = {'pD_0_/_2','pD_1_/_2','pT_0_/_2','pT^*_1_/_2','pD^*_-_/_2'};
%             elseif strcmp(prefix,'sf3')
%                 Cmap = [0 0 1; 0 0.75 0.5; 0 1 1; 1 0 0; .5 0 .75; 1 0 1];
%                 legendnames = {'pD_0_/_3','pD_1_/_3','pD_2_/_3','pT_0_/_3','pT_1_/_3','pT^*_2_/_3'};    
%             else
%                 Cmap = [0 0 1; 0 1 1; 1 0 0; 1 0 1; 0 0.75 0];
%                 legendnames = {'aD_0_/_2','aD_1_/_2','aT_0_/_2','aT^*_1_/_2','aD^*_-_/_2'};
%             end
        else
            Cmap = [1 0 0; 1 .5 0; .75 .75 0; 0 1 0; 0 0 1];
            if strcmp(prefix,'sf') || strcmp(prefix,'sf3')
                legendnames = {'pSq#1','pSq#2','pSq#3','pSq#4','pSq#5'};
            else
                legendnames = {'aSq#1','aSq#2','aSq#3','aSq#4','aSq#5'};
            end
        end

%         iFig = iFig+1;
%         figure(iFig);
%         MakeFigureTitle(sprintf('Figure %d: %s-%s%s',gcf,prefix,bandname,glmType),0);
        figure;
        PlotGroupSvdResults(avgW{iGlmType,iBand}, avgTC{iGlmType,iBand}, ...
            chanlocs{iGlmType,iBand}, tResponse{iGlmType,iBand}, ...
            steW{iGlmType,iBand}, steTC{iGlmType,iBand}, legendnames, Cmap);
        MakeFigureTitle(sprintf('Figure %d: %s-%s%s',gcf,prefix,bandname,glmType),0);
    end
end

%% Run stats
for iGlmType = 1:numel(glmTypes)
    glmType = glmTypes{iGlmType};
    for iBand = 1%:numel(bandnames);

        bandname = bandnames{iBand};
        fprintf('-----Plotting SVD stats for %s glmType %s, band %s\n',prefix, glmType,bandname);

        if ~exist(sprintf('%s-%d-GLMresults-%s%s.mat',prefix,subjects(1),bandname,glmType),'file')
            disp('SKIPPING!');
            continue;
        end
        
        if strcmp(glmType(end-9:end-6),'Type') 
            if strcmp(prefix,'sf')
                Cmap = [0 0 1; 0 1 1; 1 0 0; 1 0 1; 0 0.75 0];
                legendnames = {'pD_0_/_2','pD_1_/_2','pT_0_/_2','pT^*_1_/_2','pD^*_-_/_2'};
            elseif strcmp(prefix,'sf3')
                Cmap = [0 0 1; 0 0.75 0.5; 0 1 1; 1 0 0; .5 0 .75; 1 0 1];
                legendnames = {'pD_0_/_3','pD_1_/_3','pD_2_/_3','pT_0_/_3','pT_1_/_3','pT^*_2_/_3'};    
            else
                Cmap = [0 0 1; 0 1 1; 1 0 0; 1 0 1; 0 0.75 0];
                legendnames = {'aD_0_/_2','aD_1_/_2','aT_0_/_2','aT^*_1_/_2','aD^*_-_/_2'};
            end
        else
            Cmap = [1 0 0; 1 .5 0; .75 .75 0; 0 1 0; 0 0 1];
            if strcmp(prefix,'sf') || strcmp(prefix,'sf3')
                legendnames = {'pSq#1','pSq#2','pSq#3','pSq#4','pSq#5'};
            else
                legendnames = {'aSq#1','aSq#2','aSq#3','aSq#4','aSq#5'};
            end
        end
        
        % Get & Plot
        figure;
        GetCrossSubjectStats(timecourse{iGlmType,iBand},tResponse{iGlmType,iBand},legendnames,Cmap);
        MakeFigureTitle(sprintf('Figure %d: %s-%s%s',gcf,prefix,bandname,glmType),0);
    end
end
