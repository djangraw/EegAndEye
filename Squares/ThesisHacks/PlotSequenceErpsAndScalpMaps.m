% function PlotSequenceErpsAndScalpMaps
% Created 4/19/14 by DJ.

prefixes = {'sq','sf','sf3'};
events = {'SqNum2','SqNum1';...
    'sf-SqNum2','sf-SqNum1';...
    'sf-SqNum2','sf-SqNum1'};
event_weights = [1 -1];
baseline_win = [0 -1];
iLevel = 3;
[contrastFns, contrastVar, contrastZ] = deal(cell(1,3));

for i=1:3

    eval(sprintf('results = R_%s_sqnum;',prefixes{i}));
    event_list = events(i,:);
    
    [contrastFns{i}, contrastVar{i}, contrastZ{i}]  = SetUpTopLevelGlm_flex(results,event_list,event_weights,baseline_win,iLevel);

end

%% Run Top-Level GLMs
multcorrect = 'fdr';
[group_RF, group_P] = deal(cell(1,3));
for i=1:3
    % run level 2
    [group_RF{i},group_P{i}] = RunTopLevelGlm_EEG(contrastFns{i},contrastVar{i},multcorrect);
end

%% Get ERPs
iTimes = 1:51;
group_RF_all = [];
group_P_all = [];
for i=1:3
    group_RF_all(:,:,i) = group_RF{i}(:,iTimes);
    group_P_all(:,:,i) = group_P{i}(:,iTimes);
end
% tResponse = R_sf_sqnum(1).tResponse{end}(iTimes);
tResponse = R_sf_sqnum(1).tResponse(iTimes);
chanlocs = R_sf_sqnum(1).EEG.chanlocs;

%% Save results
save TopLevelGlmResults_Sq2vs1 group_* contrast* events event_weights baseline_win iLevel multcorrect tResponse chanlocs 

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
    ylim([-2 4]) % to match target plot
    title('');
end
%% Plot Scalp Maps
tBinCenters = 25:25:475;%[125 175 325];%37.5:75:500;
tBinWidth = 50;%75;
clim = [-2.65 2.65]; % to match target plot
figure;
% clf;
set(gcf,'Position',[684 1001 447 428]);

[sm_all,sm] = GetScalpMaps(group_RF_all,tResponse,tBinCenters,tBinWidth);
PlotScalpMaps(sm_all,chanlocs,clim,tBinCenters-tBinWidth/2,legendstr);



%% Plot Z score Scalp Maps
tBinCenters = 25:25:475;%[125 175 325];%37.5:75:500;
tBinWidth = 50;%75;
cthresh = 1.96; % z score for 2-tailed p=0.05
clim = [-5 5];
group_Z_all = norminv(group_P_all);

figure;
% clf;
set(gcf,'Position',[684 1001 447 428]);

[sm_all,sm] = GetScalpMaps(group_Z_all,tResponse,tBinCenters,tBinWidth,cthresh);
PlotScalpMaps(sm_all,chanlocs,clim,tBinCenters-tBinWidth/2,legendstr,cthresh);

