%% RECONSTRUCT FULL TRIALS AND PLOT THEM

iEvents = 4:7;
iEvents2 = 11:14;
data = newRF(:,:,iEvents)+newRF(:,:,iEvents2);%R3.responseFns{1}(:,:,iEvents);% + R.responseFns{1}(:,:,iEvents2);
tResponse = R3.tResponse{1};
chansToPlot = {'FZ';'CZ';'PZ';'OZ'};

PlotResponseFnsGrid(data,R.regressor_events{1}(iEvents),tResponse,R.EEG.chanlocs,chansToPlot,GetSquaresEventColormap(R.regressor_events{1}(iEvents)));
%%
PlotReconstructedTrial(R.responseFns{1},R.tResponse{1},R.regressor_events{1},{'pD_{0/2}','pT_{0/2}','pD_{1/2}','pT^*_{1/2}','pD_{+/2}','sf-Circle'},R.EEG.chanlocs,chansToPlot);





%% RUN GROUP GLM ON INDIVIDUAL EVENTS
multcorrect = 'none';
for i=1:size(contrastFns_sf,3)
    [group_RF_cell{i},group_P_cell{i}] = RunTopLevelGlm_EEG(contrastFns_sf(:,:,i,:),contrastVar_sf(:,:,i,:),multcorrect);
end


%%
group_RF_sf_wascell = cat(3,group_RF_cell{:});
group_P_sf_wascell = cat(3,group_P_cell{:});





%% RUN RIDGE REGRESSION ON A GIVEN EEG SESSION

EEG = R.EEG;
method = 'ridge';
lambda = 1e4;
[newRF, newTR] = RunEegGlm(EEG,R.regressor_events{R.iLevel},R.influence{R.iLevel},R.artifact_events,R.artifact_influence{R.iLevel},method,0,lambda);
