% GetComponentContrasts_script.m
%
% Created 9/19/14 by DJ.
% Updated 1/15/15 by DJ - multi-experiment.

%% Get Stats for this comparison
experiments = {'sf'};%{'sq','sf','sf3'};%
compRule = 'T1vD0';%'T2vD0';%'T1vD0';%
suffixes = {'Type-v3pt6-10fold'};
multcorrect = 'none';
% get component scalp map
[group_Z, group_RF, legendstr, titlestr] = CompileContrasts(experiments,suffixes,{compRule},multcorrect);

% Get chanlocs and times
foo = load(sprintf('%s-%sResults-%scontrast.mat',experiments{1},suffixes{1},compRule),'chanlocs','tResponse');
chanlocs = foo.chanlocs; tResponse = foo.tResponse;

figure(11); clf;
% tLims = [100 600]; % in ms
tLims = [100 350]; % in ms
tComp = nan(1,numel(experiments));
compSM = nan(numel(chanlocs),numel(experiments));
for iExp = 1:numel(experiments);
    % get component of max power in T2vD0 contrast
    iOkTime = find(tResponse>tLims(1) & tResponse<tLims(2));
    [~,iMax] = max(sum(group_RF(:,iOkTime,iExp).^2,1));
    compSM(:,iExp) = group_RF(:,iOkTime(iMax),iExp);
    % normalize
    compSM(:,iExp) = compSM(:,iExp)/sqrt(sum(compSM(:,iExp).^2));


    subplot(1,numel(experiments),iExp);
    topoplot(compSM(:,iExp),chanlocs);
    colorbar;
    tComp(iExp) = tResponse(iOkTime(iMax));
    title(sprintf('%s Scalp map at max power time (t=%dms)',experiments{iExp},tComp(iExp)))
end
%% get component betas

foo = load('Type-v3pt6-Peak-10fold_vs_Type-v3pt6-10fold-Betas','betas','lambda_best');
exp_all = {'sq','sf','sf3'};
iExp = nan(1,numel(experiments));
for i=1:numel(experiments)
    iExp(i) = find(strcmp(experiments{i},exp_all));
end
betas = foo.betas(2,iExp);
        
lambda_best = foo.lambda_best(2,iExp);

%% Calculate component contrasts
iExp = 1; % in experiments list
old_suffix = 'GLMresults-Type-v3pt6-RampUp';
suffix_in = 'Type-v3pt6-Matrices';
suffix_out = sprintf('Type-v3pt6-%s%dmscomp',compRule,tComp(iExp));

switch experiment
    case 'sf3'
        rules = {'D2vD0','D1vD0','T0vD0','T1vD0','T2vD0',  'D2vD1','T1vT0','T2vT1'};
    otherwise
        rules = {'D1vD0','T0vD0','T1vD0','T1vT0'};
end

% ---Added 1/15/15: START
iEvents_cell = cell(1,numel(experiments));
for i=1:numel(experiments)  
    experiment = experiments{i};

        if strcmp(experiment,'sf3')
            iEvents_cell{i} = 2:14; % for basic model
        else
            iEvents_cell{i} = 2:12; % for basic model
        end        
end
% iEvents = iEvents_cell{iExp};
% ---Added 1/15/15: END
GetComponentContrast(rules,experiments,old_suffix,suffix_in,suffix_out,iEvents_cell,lambda_best,betas,compSM')

%%
multcorrect = 'none';

switch experiment
    case 'sf3'
        statrules = {'D2vD1','D1vD0','T0vD0','T1vT0','T2vT1'};
        barrules = {'D2vD0','D1vD0','T0vD0','T1vD0','T2vD0'};
    otherwise
        statrules = {'D1vD0','T0vD0','T1vT0'};
        barrules = {'D1vD0','T0vD0','T1vD0'};
end

[comp_Z, comp_RF, comp_str] = CompileContrasts(experiments,{suffix_out},statrules,multcorrect);
[bar_Z0, bar_RF, bar_str0] = CompileContrasts(experiments,{suffix_out},barrules,multcorrect);
%%
% iMax = 50;
iWin = (iOkTime(iMax)-2):(iOkTime(iMax)+2);
% [~,iMax] = max(bar_Z(1,:,end),[],2);
% iWin = (iMax-2):(iMax+2);

switch experiment
    case 'sf3'
        barevents = {'pD_{2/3}','pD_{1/3}','pD_{0/3}','pT_{0/3}','pT_{1/3}','pT^*_{2/3}'};
        iZero = find(strcmp('pD_{0/3}',barevents));
        eventsign = [-1 -1 -1 1 1 1];
    case 'sf'
        barevents = {'pD_{1/2}','pD_{0/2}','pT_{0/2}','pT^*_{1/2}'};
        iZero = find(strcmp('pD_{0/2}',barevents));
        eventsign = [-1 -1 1 1];
    case 'sq'
        barevents = {'aD_{1/2}','aD_{0/2}','aT_{0/2}','aT^*_{1/2}'};
        iZero = find(strcmp('aD_{0/2}',barevents));
        eventsign = [-1 -1 1 1];
end
% insert zeros        
bar_Z = zeros(size(bar_Z0,1),size(bar_Z0,2),size(bar_Z0,3)+1);
bar_Z(:,:,1:(iZero-1)) = bar_Z0(:,:,1:(iZero-1));
bar_Z(:,:,(iZero+1):end) = bar_Z0(:,:,iZero:end);
% flip signs on distractor comparisons
bar_Z = bar_Z .* repmat(permute(eventsign,[3 1 2]),size(bar_Z,1),size(bar_Z,2));
comp_Z = comp_Z .* repmat(permute(eventsign([1:(iZero-1), (iZero+1):end]),[3 1 2]),size(bar_Z,1),size(bar_Z,2));

bar_str = [bar_str0(1:(iZero-1)), {'zero'}, bar_str0(iZero:end)];

Cmap = GetSquaresEventColormap(barevents);

% PLOT!
clf;
subplot(3,1,1);cla;
% PlotResponseFnsGrid(bar_Z, bar_str,tResponse,chanlocs(1),{chanlocs(1).labels},Cmap);
PlotResponseFns(bar_Z,bar_str,tResponse,chanlocs(1),1,Cmap)
PlotVerticalLines(tResponse(iWin([1 end])),'k:');
title('Component Activity relative to D0');

subplot(3,1,2);cla;
% PlotResponseFnsGrid(comp_Z, comp_str,tResponse,chanlocs(1),{chanlocs(1).labels},Cmap);
PlotResponseFns(comp_Z,comp_str,tResponse,chanlocs(1),1,Cmap)
PlotVerticalLines(tResponse(iWin([1 end])),'k:');
title('Component Activity relative to previous event');

subplot(3,1,3);cla;
cla; hold on;
for i=1:size(bar_Z,3)
%     % use Stouffer's method to average across a window   
%     bar_newZ = sum(bar_Z(1,iWin,i),2)/sqrt(numel(iWin));    
%     bar(i,bar_newZ,'facecolor',Cmap(i,:));
    bar(i,squeeze(mean(bar_Z(:,iWin,i),2)),'facecolor',Cmap(i,:));
    
    if i<size(bar_Z,3)
%         % Use Stouffer's method
%         comp_newZ = sum(comp_Z(1,iWin,i),2)/sqrt(numel(iWin));
%         plot(i+.5,comp_newZ,'*','color',Cmap(i,:));
        plot(i+.5,squeeze(mean(comp_Z(:,iWin,i),2)),'*','color',Cmap(i,:));
    end    
end
PlotHorizontalLines([-1.97 1.97],'k--');
set(gca,'xtick',1:size(bar_Z,3),'xticklabel',barevents);
xlabel('Event');
ylabel('Mean Z score in window');
title('Z scores in window around max');
MakeFigureTitle(sprintf('%s experiment',experiment));
