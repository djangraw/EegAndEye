%% Set up
[RF_cell,CM_cell] = deal(cell(1,numel(R_sf3_type)));
T = size(R_sf3_type(1).responseFns{3},2);

%% Set up Top-Level SF3 GLM
% event_list_sf3 = {'pT_{0/3}','pT_{1/3}','pT^*_{2/3}'};
% event_weights_sf3 = [1 1 1;... % is a target
%                     1/3 2/3 1;... % level of evidence after this stimulus
%                     0 0 1]; % is a decision
event_list_sf3 = {'pD_{0/3}','pT_{0/3}','pT_{1/3}','pT^*_{2/3}'};
event_weights_sf3 = [-1 1 0 0;... % T0-D0
                    -1 0 1 0;... % T1-D0
                    -1 0 0 1]; % T2-D0              
% event_list_sf3 = {'pD_{0/3}','pT_{0/3}','pT_{1/3}','pT^*_{2/3}'};
% event_weights_sf3 = [-1 1/3 1/3 1/3;... % Target-D0
%                     -1 1/6 2/6 3/6;... % evidence-D0
%                     -1 0 0 1]; % decision-D0     
% % normalize
% for i=1:size(event_weights_sf3,1);
%     event_weights_sf3(i,:) = event_weights_sf3(i,:)/sum(abs(event_weights_sf3(i,:)))/2;
% end
[contrastFns_sf3, contrastVar_sf3, contrastZ_sf3]  = SetUpTopLevelGlm_flex(R_sf3_type,event_list_sf3,event_weights_sf3);

%% Run Top-Level SF3 GLM
multcorrect = 'none';
% run level 2
[group_RF_sf3,group_P_sf3] = RunTopLevelGlm_EEG(contrastFns_sf3,contrastVar_sf3,multcorrect);



%% Set up Top-level SF GLM
% event_list_sf = {'pT_{0/2}','pT^*_{1/2}'};
% event_weights_sf = [1 1;... % is a target
%                     0 1]; % is a decision
event_list_sf = {'pD_{0/2}','pT_{0/2}','pT^*_{1/2}'};
event_weights_sf = [-1 1 0;... % T0-D0
                    -1 0 1]; % T2-D0
% event_list_sf = {'pD_{0/2}','pT_{0/2}','pT^*_{1/2}'};
% event_weights_sf = [-1 1/2 1/2;... % target-D0
%                     -2 2/5 8/5]; % evidence-D0 + decision-D0                
% %normalize
% for i=1:size(event_weights_sf,1);
%     event_weights_sf(i,:) = event_weights_sf(i,:)/sum(abs(event_weights_sf(i,:)))/2;
% end
% run level 2
[contrastFns_sf, contrastVar_sf, contrastZ_sf]  = SetUpTopLevelGlm_flex(R_sf_type,event_list_sf,event_weights_sf);

%% Run Top-level SF GLM
multcorrect = 'none';
[group_RF_sf,group_P_sf] = RunTopLevelGlm_EEG(contrastFns_sf,contrastVar_sf,multcorrect);



%% Set up Top-level SQ GLM
% event_list_sq = {'aT_{0/2}' 'aT^*_{1/2}'};
% event_weights_sq = [1 1;... % is a target
%                     0 1]; % is a decision
event_list_sq = {'aD_{0/2}','aT_{0/2}','aT^*_{1/2}'};
event_weights_sq = [-1 1 0;... % is a target
                    -1 0 1]; % is a decision
% event_list_sq = {'aD_{0/2}','aT_{0/2}','aT^*_{1/2}'};
% event_weights_sq = [-1 1/2 1/2;... % target-D0
%                     -2 2/5 8/5]; % evidence-D0 + decision-D0    

% %normalize
% for i=1:size(event_weights_sq,1);
%     event_weights_sq(i,:) = event_weights_sq(i,:)/sum(abs(event_weights_sq(i,:)));
% end
% run level 2
[contrastFns_sq, contrastVar_sq, contrastZ_sq]  = SetUpTopLevelGlm_flex(R_sq_type,event_list_sq,event_weights_sq);

%% Run Top-level SQ GLM
multcorrect = 'none';
[group_RF_sq,group_P_sq] = RunTopLevelGlm_EEG(contrastFns_sq,contrastVar_sq,multcorrect);



%% Plot scalp maps
tBinCenters = 37.5:75:750;
tBinWidth = 75;
legendstr_sf3 = {'Target','Evidence','Decision'};
legendstr_sf = {'Target','Evidence + Decision'};
legendstr_sq = {'Target','Evidence + Decision'};

figure(1); clf;
[sm_all,sm] = GetScalpMaps(group_RF_sf3,R_sf3_type(1).tResponse{end},tBinCenters,tBinWidth);
PlotScalpMaps(sm_all,R_sf3_type(1).EEG.chanlocs,[],tBinCenters-tBinWidth/2,legendstr_sf3);
% Annotate figure
MakeFigureTitle(sprintf('Level 3 results for %d sf3 Subjects',numel(R_sf3_type)));

figure(2); clf;
[sm_all,sm] = GetScalpMaps(group_RF_sf,R_sf_type(1).tResponse{end},tBinCenters,tBinWidth);
PlotScalpMaps(sm_all,R_sf_type(1).EEG.chanlocs,[],tBinCenters-tBinWidth/2,legendstr_sf);
% Annotate figure
MakeFigureTitle(sprintf('Level 3 results for %d sf Subjects',numel(R_sf_type)));

figure(5); clf;
[sm_all,sm] = GetScalpMaps(group_RF_sq,R_sq_type(1).tResponse{end},tBinCenters,tBinWidth);
PlotScalpMaps(sm_all,R_sq_type(1).EEG.chanlocs,[],tBinCenters-tBinWidth/2,legendstr_sq);
% Annotate figure
MakeFigureTitle(sprintf('Level 3 results for %d sq Subjects',numel(R_sq_type)));


%% Plot statistics
cthresh = 1.96; % z score for 2-tailed p=0.05

figure(3); clf;
group_Z_sf3 = norminv(group_P_sf3);
[sm_all,sm] = GetScalpMaps(group_Z_sf3,R_sf3_type(1).tResponse{end},tBinCenters,tBinWidth,cthresh);
PlotScalpMaps(sm_all,R_sf3_type(1).EEG.chanlocs,[],tBinCenters-tBinWidth/2,legendstr_sf3);
% Annotate figure
MakeFigureTitle(sprintf('Level 3 results for %d sf3 Subjects',numel(R_sf3_type)));
%%
figure(4); clf;
group_Z_sf = norminv(group_P_sf);
[sm_all,sm] = GetScalpMaps(group_Z_sf,R_sf_type(1).tResponse{end},tBinCenters,tBinWidth,cthresh);
PlotScalpMaps(sm_all,R_sf_type(1).EEG.chanlocs,[],tBinCenters-tBinWidth/2,legendstr_sf);
% Annotate figure
MakeFigureTitle(sprintf('Level 3 results for %d sf Subjects',numel(R_sf_type)));

figure(6); clf;
group_Z_sq = norminv(group_P_sq);
[sm_all,sm] = GetScalpMaps(group_Z_sq,R_sq_type(1).tResponse{end},tBinCenters,tBinWidth,cthresh);
PlotScalpMaps(sm_all,R_sq_type(1).EEG.chanlocs,[],tBinCenters-tBinWidth/2,legendstr_sq);
% Annotate figure
MakeFigureTitle(sprintf('Level 3 results for %d sq Subjects',numel(R_sq_type)));


%% TRY GROUP X MATRIX AND CONTRASTS!!!
% sf3: target, evidence, decision
[D,T,M,N] = size(contrastFns_sf3);
contrastFns_all = reshape(contrastFns_sf3,[D T M*N]);
contrastVar_all = reshape(contrastVar_sf3,[D T M*N]);
% group_x = repmat([1 0 0; 0 1 0; 0 0 1],N,1);
group_x = repmat([1 0 0; 1 1/2 0; 1 1 1],N,1); % Rows = T0,T1,T2; Cols = T,E,D

% sf: target, evidence + decision
[D,T,M,N] = size(contrastFns_sf);
contrastFns_all = cat(3,contrastFns_all,reshape(contrastFns_sf,[D T M*N]));
contrastVar_all = cat(3,contrastVar_all,reshape(contrastVar_sf,[D T M*N]));
% group_x = cat(1,group_x,repmat([1 0 0; 0 1 1],N,1)); 
group_x = cat(1,group_x,repmat([1 0 0; 1 1 1],N,1)); % Rows = T0,T1; Cols = T,E,D

% sq: target, evidence + decision
[D,T,M,N] = size(contrastFns_sq);
contrastFns_all = cat(3,contrastFns_all,reshape(contrastFns_sq,[D T M*N]));
contrastVar_all = cat(3,contrastVar_all,reshape(contrastVar_sq,[D T M*N]));
group_x = cat(1,group_x,repmat([1 0 0; 0 1 1],N,1));

group_contrast = eye(size(group_x,2));
%% Run top-level GLM with group x and contrasts
multcorrect = 'none';
[group_RF, group_p] = RunTopLevelGlm_EEG_flex(contrastFns_all,contrastVar_all,group_x,group_contrast,multcorrect);

%% Plot results

%% Plot scalp maps
tBinCenters = 37.5:75:750;
tBinWidth = 75;
cthresh = 1.96; % z score for 2-tailed p=0.05
legendstr_all = {'Target','Evidence','Decision'};
chanlocs = R_sf3_type(1).EEG.chanlocs;
tResponse = R_sf3_type(1).tResponse{end};

figure(9); clf;
[sm_all,sm] = GetScalpMaps(group_RF,tResponse,tBinCenters,tBinWidth);
PlotScalpMaps(sm_all,chanlocs,[],tBinCenters-tBinWidth/2,legendstr_all);
% Annotate figure
MakeFigureTitle(sprintf('Level 3 results for all experiments (%d Subjects)',numel(R_sf3_type)+numel(R_sf_type)+numel(R_sq_type)));
% MakeFigureTitle(sprintf('Level 3 results for all passive experiments (%d Subjects)',numel(R_sf3_type)+numel(R_sf_type)));
% MakeFigureTitle(sprintf('Level 3 results for all SF3 experiments (%d Subjects)',numel(R_sf3_type)));
figure(10); clf;
group_Z = norminv(group_p);
[sm_all,sm] = GetScalpMaps(group_Z,tResponse,tBinCenters,tBinWidth,cthresh);
PlotScalpMaps(sm_all,chanlocs,[],tBinCenters-tBinWidth/2,legendstr_all);
% Annotate figure
MakeFigureTitle(sprintf('Level 3 results for all experiments (%d Subjects)',numel(R_sf3_type)+numel(R_sf_type)+numel(R_sq_type)));
% MakeFigureTitle(sprintf('Level 3 results for all passive experiments (%d Subjects)',numel(R_sf3_type)+numel(R_sf_type)));
% MakeFigureTitle(sprintf('Level 3 results for all SF3 experiments (%d Subjects)',numel(R_sf3_type)));