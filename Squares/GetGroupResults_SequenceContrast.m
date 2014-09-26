% GetGroupResults_SequenceContrast
%
% This script builds a sequence-based contrast matrix, calculates the
% individual and group RFs and stats for the contrast, and plots them.
%
% Created 8/11/14 by DJ.

%% Load data
% experiment = 'sf';
% if strcmp(experiment,'sf')
%     subjects = [1:10 12:13];
% elseif strcmp(experiment,'sq')
%     subjects = [9:11 13:15 17:19 20:27];
% else
%     subjects = 1:12;
% end
% % suffix = '-Type-v3pt6-RampUp-ridge100';
% suffix = '-Type-v3pt6-RampUp-ridge014';
% % suffix = '-v3pt5-RampUp';
% 
% clear R
% for iSubj=1:numel(subjects)
%     % Set up
%     fprintf('Loading %s-%d-GLMresults%s...\n',experiment,subjects(iSubj),suffix);
%     R(iSubj) = load(sprintf('%s-%d-GLMresults%s',experiment,subjects(iSubj),suffix));
% end


%%

event_list = R(1).regressor_events{end};
tResponse = R(1).tResponse{end};
% rule = 'T0vD0';
if strcmp(experiment,'sq')
    tEvents = 0:1000:3750;
else
    tEvents = 0:500:2500;
end


if strcmp(experiment,'sf3')
    switch rule
        case 'T1vT0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pT_{0/3}','pD^*_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pD^*_{-/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdTdd-ddTdd'};
        case 'T2vT1'
            sequence{1} = {'pD_{0/3}','pT_{0/3}','pT_{1/3}','pD_{2/3}','pD^*_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pD_{+/3}','pD_{+/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'ttTdd-dtTdd'};
        case 'T2vT0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pT_{0/3}','pD^*_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pD_{+/3}','pD_{+/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'ttTdd-ddTdd'};            
            
        case 'T0vD0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pD_{0/3}','pD_{0/3}','pT_{0/3}','pD^*_{-/3}','pD_{-/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'ddTdd-dDddd'};
        case 'D1vD0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pD_{1/3}','pD_{1/3}','pD_{1/3}','pD^*_{-/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdDdd-dDddd'};
        case 'T1vD0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pD^*_{-/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdTdd-dDddd'};
        case 'D2vD0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pT_{1/3}','pD_{2/3}','pD_{2/3}','pD^*_{-/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(4);
            legendstr = {'ttDdd-dDddd'}; 
        case 'T2vD0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pD_{+/3}','pD_{+/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'ttTdd-dDddd'}; 
        otherwise
            sequence = {{},{}};
    end
else


    switch rule
        case 'TDvDD'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','NULL','sf-Circle'};
            sequence{2} = {'pT_{0/2}','pD_{1/2}','pD_{1/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};        
            tContrast{1} = tResponse + tEvents(1);
            tContrast{2} = tResponse + tEvents(1);
            legendstr = {'Tdddd-Ddddd'};
        case 'T1vT0'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pT_{0/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};
            sequence{2} = {'pT_{0/2}','pD_{1/2}','pT^*_{1/2}','pD_{+/2}','pD_{+/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdTdd-ddTdd'};
        case 'T1vD0'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','NULL','sf-Circle'};
            sequence{2} = {'pT_{0/2}','pD_{1/2}','pT^*_{1/2}','pD_{+/2}','pD_{+/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdTdd-ddDdd'};
        case 'T0sn4vT0sn2'
            sequence{1} = {'pD_{0/2}','pT_{0/2}','pD_{1/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};
            sequence{2} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pT_{0/2}','pD^*_{-/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(4);
            legendstr = {'dddTd-dTddd'};
        case 'T1sn4vT1sn2'
            sequence{1} = {'pT_{0/2}','pT^*_{1/2}','pD_{+/2}','pD_{+/2}','pD_{+/2}','sf-Circle'};
            sequence{2} = {'pD_{0/2}','pD_{0/2}','pT_{0/2}','pT^*_{1/2}','pD_{+/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(4);
            legendstr = {'ddtTd-tTddd'};
        case 'D1vD0'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','NULL','sf-Circle'};
            sequence{2} = {'pT_{0/2}','pD_{1/2}','pD_{1/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdDdd-ddDdd'};
        case 'T0vD0'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','NULL','sf-Circle'};
            sequence{2} = {'pD_{0/2}','pD_{0/2}','pT_{0/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'ddTdd-ddDdd'};
        otherwise
            sequence = {{},{}};
    end
end

if strcmp(experiment,'sq')
    for i=1:2
        for j=1:numel(sequence{i})
            if sequence{i}{j}(1) == 'p'
                sequence{i}{j}(1) = 'a';
            end
        end
    end
end

C_minus = BuildSequenceContrast(event_list,tResponse,sequence{1},tEvents,tContrast{1});
C_plus = BuildSequenceContrast(event_list,tResponse,sequence{2},tEvents,tContrast{2});

contrasts = C_plus-C_minus;
% Plot contrast matrix
figure(55);
imagesc(tResponse,1:size(contrasts,1),contrasts);
xlabel('tResponse')
ylabel('Event')
T = length(tResponse);
set(gca,'ytick',(1:T:size(contrasts,1)),'yticklabel',event_list);
%% Get single-subject contrasts
% set up contrast functions
N = numel(R);
D = size(R(1).responseFns{R(1).iLevel},1);
T = size(R(1).responseFns{R(1).iLevel},2);
M = size(contrasts,2)/T; % number of contrasts
% [contrastFns, contrastVar, contrastZ] = deal(nan(D,T,M,N));
doPlot = false;
multcorrect = 'none';
% find contrast functions
for i=1:N
    fprintf('---subject %d/%d...\n',i,N);
    % prepare figure    
    if doPlot
        figure(100+i); clf;  
        MakeFigureTitle(sprintf('%s, unknown contrast',R(i).filenames{R(i).iLevel}),1);
    end
    % Calculate and plot
    [contrastZ(:,:,:,i),~,~,contrastFns(:,:,:,i),contrastVar(:,:,:,i)] = GetGlmZscore(R(i).EEG,R(i),multcorrect,contrasts,doPlot);
end

%% Subtract baseline
% contrastFns0 = contrastFns;
% iBaseline = 1:10;
% contrastFns = contrastFns0 - repmat(mean(contrastFns0(:,iBaseline,:,:),2),[1,101,1,1]);


%% Plot responses
% subjstr = cell(1,N);
% for i=1:N
%     subjstr{i} = sprintf('Subj %d',subjects(i));
% end
% figure(201); clf;
% [sm_all,sm] = GetScalpMaps(permute(contrastFns,[1 2 4 3]),tResponse,tBinCenters,tBinWidth);
% PlotScalpMaps(sm_all,chanlocs,[],tBinCenters-tBinWidth/2,subjstr);
% MakeFigureTitle(sprintf('%s %s, %s contrast, Single-Subject RFs',experiment,suffix,rule));
% 
% %% Plot Z scores
% cthresh = 1.96; % z score for 2-tailed p=0.05
% figure(202); clf;
% [sm_all,sm] = GetScalpMaps(permute(contrastFns,[1 2 4 3]),tResponse,tBinCenters,tBinWidth,cthresh);
% PlotScalpMaps(sm_all,chanlocs,[],tBinCenters-tBinWidth/2,subjstr);
% MakeFigureTitle(sprintf('%s %s, %s contrast, Single-Subject Z scores',experiment,suffix,rule));


%% Get group-level stats
multcorrect = 'none';
[group_RF,group_P] = RunTopLevelGlm_EEG(contrastFns,contrastVar,multcorrect);
group_Z = norminv(group_P);

%% Plot group results as ERPs


% chansToPlot = {'FZ';'CZ';'PZ';'OZ'};
% % Cmap = GetSquaresEventColormap({'pT_{0/2}'});
% % Cmap = GetSquaresEventColormap({'pT_{0/2}','pD_{1/2}','pT^*_{1/2}'});
% Cmap = GetSquaresEventColormap({'pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pT^*_{2/3}'});
% % Cmap = GetSquaresEventColormap(event_list_sf(2:6));
% 
% figure(61); clf;
% PlotResponseFnsGrid(group_RF, legendstr,tResponse,chanlocs,chansToPlot,Cmap);
% % MakeFigureTitle(suff);
% set(gcf,'Position',[1 268 568 1238]);
% for i=1:4
% subplot(4,1,i);
% title('')
% ylabel(chansToPlot{i});
% end
% subplot(4,1,1);
% title(sprintf('%s %s, Group RFs',experiment,suffix))
% 
% figure(71); clf;
% PlotResponseFnsGrid(group_Z, legendstr,tResponse,chanlocs,chansToPlot,Cmap);
% % MakeFigureTitle(suff);
% set(gcf,'Position',[569 268 568 1238]);
% for i=1:4
% subplot(4,1,i);
% title('')
% ylabel(chansToPlot{i});
% end
% subplot(4,1,1);
% title(sprintf('%s %s, Group Z scores',experiment,suffix))
% 
% %% Plot scalp maps
% tBinCenters = 37.5:75:750;
% tBinWidth = 75;
% 
% figure(2); clf;
% [sm_all,sm] = GetScalpMaps(group_RF,tResponse,tBinCenters,tBinWidth);
% PlotScalpMaps(sm_all,chanlocs,[],tBinCenters-tBinWidth/2,legendstr);
% % Annotate figure
% MakeFigureTitle(sprintf('Group RFs for %d %s Subjects: %s',size(contrastFns,4),experiment,legendstr{1}));
% % set(gcf,'Position',[825 1352 1261 145]);
% 
% 
% % Plot statistics
% cthresh = 1.96; % z score for 2-tailed p=0.05
% 
% figure(3); clf;
% % group_Z = norminv(group_P);
% [sm_all,sm] = GetScalpMaps(group_Z,tResponse,tBinCenters,tBinWidth,cthresh);
% PlotScalpMaps(sm_all,chanlocs,[],tBinCenters-tBinWidth/2,legendstr);
% % Annotate figure
% MakeFigureTitle(sprintf('Group Z scores for %d %s Subjects: %s',size(contrastFns,4),experiment,legendstr{1}));
% % set(gcf,'Position',[825 1202 1261 145]);


%% Save results
save(sprintf('%s%sResults-%scontrast',experiment,suffix,rule),'event_*','contrast*',...
    'group_*','chanlocs','tResponse','experiment','subjects','suffix',...
    'rule','sequence','tContrast','legendstr')

%% Reconstruct multiple contrasts
% if strcmp(experiment,'sf3')
%     rules = {'T0vD0','D1vD0','T1vD0','D2vD0','T2vD0'};
% else
%     rules = {'T0vD0','D1vD0','T1vD0'};
% end
% [group_RF,group_Z,contrastFns] = deal([]);
% legendstr = {};
% for i=1:numel(rules)
%     foo = load(sprintf('%s%sResults-%scontrast',experiment,suffix,rules{i}));
%     group_RF = cat(3,group_RF,foo.group_RF);
%     group_Z = cat(3,group_Z,foo.group_Z);
%     contrastFns = cat(3,contrastFns,foo.contrastFns);
%     legendstr = cat(2,legendstr,foo.legendstr);
% end