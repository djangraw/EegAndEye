%% Get Stats for this comparison
experiments = {'sf3'};
compRule = 'T2vD0';
suffixes = {'Type-v3pt6-10fold'};
% get component scalp map
[group_Z, group_RF, legendstr, titlestr] = CompileContrasts(experiments,suffixes,{compRule},multcorrect);

% get component of max power in T2vD0 contrast
[~,iMax] = max(sum(group_RF(:,1:40).^2,1));
compSM = group_RF(:,iMax);
% normalize
compSM = compSM/sqrt(sum(compSM.^2));

figure;
topoplot(compSM,chanlocs);
colorbar;
title(sprintf('Scalp map at max power time (t=%dms)',tResponse(iMax)))
%% get component betas

foo = load('Type-v3pt6-Peak-10fold_vs_Type-v3pt6-10fold-Betas','betas','lambda_best');
iExp = 3;

betas = foo.betas(2,iExp);
        
lambda_best = foo.lambda_best(2,iExp);

%% Calculate component contrasts
old_suffix = 'GLMresults-Type-v3pt6-RampUp';
suffix_in = 'Type-v3pt6-Matrices';
suffix_out = 'Type-v3pt6-T2vD0270mscomp';


rules = {'D2vD0','D1vD0','T0vD0','T1vD0','T2vD0',  'D2vD1','T1vT0','T2vT1'};

GetComponentContrast(rules,experiments,old_suffix,suffix_in,suffix_out,iEvents,lambda_best,betas,compSM')

%%
multcorrect = 'none';

statrules = {'D2vD1','D1vD0','T0vD0','T1vT0','T2vT1'};
[comp_Z, comp_RF, comp_str] = CompileContrasts(experiments,{suffix_out},statrules,multcorrect);

barrules = {'D2vD0','D1vD0','T0vD0','T1vD0','T2vD0'};
[bar_Z0, bar_RF, bar_str] = CompileContrasts(experiments,{suffix_out},barrules,multcorrect);
%%
% iMax = 50;
iWin = (iMax-2):(iMax+2);
% [~,iMax] = max(bar_Z(1,:,end),[],2);
% iWin = (iMax-2):(iMax+2);

% insert zeros
bar_Z = zeros(size(bar_Z0,1),size(bar_Z0,2),size(bar_Z0,3)+1);
bar_Z(:,:,1:2) = bar_Z0(:,:,1:2);
bar_Z(:,:,4:end) = bar_Z0(:,:,3:end);

barevents = {'pD_{2/3}','pD_{1/3}','pD_{0/3}','pT_{0/3}','pT_{1/3}','pT^*_{2/3}'};
Cmap = GetSquaresEventColormap(barevents);


subplot(3,1,1);cla;
% PlotResponseFnsGrid(bar_Z, bar_str,tResponse,chanlocs(1),{chanlocs(1).labels},Cmap);
PlotResponseFns(bar_Z,bar_str,tResponse,chanlocs(1),1,Cmap)

subplot(3,1,2);cla;
% PlotResponseFnsGrid(comp_Z, comp_str,tResponse,chanlocs(1),{chanlocs(1).labels},Cmap);
PlotResponseFns(comp_Z,comp_str,tResponse,chanlocs(1),1,Cmap)

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
set(gca,'xtick',1:size(bar_Z,3),'xticklabel',barevents);

