% CompareResidualsAcrossAnalyses
%
% Created 7/31/14 by DJ.
experiment = 'sf';
subjects = [1:10 12:13];
suffixes = {'-RampUp','-RampDown','-Peak','-Valley'};

for iSuffix = 1:numel(suffixes)
    for iSubj=1:numel(subjects)
        % Set up
        fprintf('Loading %s-%d-GLMresults-Events-Type-v3pt4%s...\n',experiment,subjects(iSubj),suffixes{iSuffix});
        R = load(sprintf('%s-%d-GLMresults-Events-Type-v3pt4%s',experiment,subjects(iSubj),suffixes{iSuffix}));
        % Subtract out
        NewEEG = R.EEG;
        for i=1:R.iLevel
            [NewEEG,S] = SubtractOutGlmResponses(NewEEG,R.responseFns{i},R.influence{i},R.regressor_events{i},R.stddev);
        end

        % Get residuals (from GetGlmZscore.m)
        Eres = reshape(NewEEG.data,NewEEG.nbchan,NewEEG.pnts*NewEEG.trials);
        isInPlay = any(S~=0,2);
        X = double(S(isInPlay,:)); % cropped regressor matrix
        [n,p] = size(X);
        fprintf('Calcluating inverse...');
        tic;
        X_inv = (X'*X)^(-1);
        t = toc;
        fprintf(' Done (took %.1f seconds).\n',t);
        r = Eres(:,isInPlay)'; % residuals
        sigmasq = sum(r.^2); % sum of squared residuals
        dof = n-p; % degrees of freedom
        % if isempty(contrasts)
            contrasts = eye(p);
        % end

        % correction_fact = contrasts'*X_inv*contrasts; % correction factor
        % beta_var = (sigmasq/dof)*correction_fact;
        % residStdErr = sqrt(beta_var);

        Nh = size(R.responseFns{R.iLevel},2);
        Nc = size(contrasts,2)/Nh; % # contrasts
        D = size(sigmasq,2); % # electrodes
        beta = reshape(R.responseFns{R.iLevel},[D p]);

        correction_fact = zeros(1,Nc*Nh);
        contrastFns = zeros(D,Nh,Nc); % like responseFns
        contrastVar = zeros(D,Nh,Nc);
        for i=1:(Nc*Nh)
            iCF = ceil(i/Nh); % which contrast fn?
            iT = i-Nh*(iCF-1); % which time point in that CF?
            correction_fact(i) = contrasts(:,i)'*X_inv*contrasts(:,i); % correction factor to account for linear dependence between contrasts
        %     contrastVar(:,iT,iCF) = sigmasq'*correction_fact(i); % variance of contrast function
            contrastVar(:,iT,iCF) = sigmasq'*correction_fact(i)/dof; % variance of contrast function
            % Get response fn of contrast    
            contrastFns(:,iT,iCF) = beta*contrasts(:,i); % convert response functions to contrast functions
        end
        % residStdErr = sqrt(contrastVar/dof); % take sqrt to get std error (???) of residuals
        residStdErr = sqrt(contrastVar); % take sqrt to get std error (???) of residuals

        subjStdErr(iSuffix,iSubj) = mean(residStdErr(:));

        fprintf('=== Mean RMS Error is %g.\n',mean(residStdErr(:)));
    end
end

%% Print results
for iSuffix = 1:numel(suffixes)
    fprintf('===== %s: Mean StdErr across %d subjects is %g.\n',suffixes{iSuffix}, numel(subjects),mean(subjStdErr(iSuffix,:)))
end

%% Load results
% suffixes = {'-v3pt5-RampUp'};
% suffixes = {'-v3pt6-RampUp-ridge100'};

for iSuffix = 1:numel(suffixes)

    for iSubj=1:numel(subjects)
        % Set up
%         fprintf('Loading %s-%d-GLMresults-Events-Type-v3pt4%s...\n',experiment,subjects(iSubj),suffixes{iSuffix});
%         R(iSubj,iSuffix) = load(sprintf('%s-%d-GLMresults-Events-Type-v3pt4%s',experiment,subjects(iSubj),suffixes{iSuffix}));
        fprintf('Loading %s-%d-GLMresults-Type%s...\n',experiment,subjects(iSubj),suffixes{iSuffix});
        R(iSubj,iSuffix) = load(sprintf('%s-%d-GLMresults-Type%s',experiment,subjects(iSubj),suffixes{iSuffix}));
    end
end

%% Get Stats
iSuffix = 1;
suff = suffixes{iSuffix};
event_list_sf = {'pD_{0/2}','pT_{0/2}','pD_{1/2}','pT^*_{1/2}',sprintf('pD_{0/2}%s',suff),sprintf('pT_{0/2}%s',suff),sprintf('pD_{1/2}%s',suff),sprintf('pT^*_{1/2}%s',suff)};
% event_list_sf = {'aD_{0/2}','aT_{0/2}','aD_{1/2}','aT^*_{1/2}',sprintf('aD_{0/2}%s',suff),sprintf('aT_{0/2}%s',suff),sprintf('aD_{1/2}%s',suff),sprintf('aT^*_{1/2}%s',suff)};
SqNum = 3;
event_weights_sf = [-1 1 0 0 -SqNum/3 SqNum/3 0 0;
                    -1 0 1 0 -SqNum/3 0 SqNum/3 0;
                    -1 0 0 1 -SqNum/3 0 0 SqNum/3];

% event_list_sf = {'pD_{0/3}','pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pT^*_{2/3}',sprintf('pD_{0/3}%s',suff),sprintf('pT_{0/3}%s',suff),sprintf('pD_{1/3}%s',suff),sprintf('pT_{1/3}%s',suff),sprintf('pD_{2/3}%s',suff),sprintf('pT^*_{2/3}%s',suff)};                
% event_weights_sf = [-1 1 0 0 0 0 -SqNum/3 SqNum/3 0 0 0 0;
%                     -1 0 1 0 0 0 -SqNum/3 0 SqNum/3 0 0 0;
%                     -1 0 0 1 0 0 -SqNum/3 0 0 SqNum/3 0 0;
%                     -1 0 0 0 1 0 -SqNum/3 0 0 0 SqNum/3 0;
%                     -1 0 0 0 0 1 -SqNum/3 0 0 0 0 SqNum/3];
                
% event_list_sf = {sprintf('pD_{0/2}%s',suff),sprintf('pT_{0/2}%s',suff),sprintf('pD_{1/2}%s',suff),sprintf('pT^*_{1/2}%s',suff)};

% event_list_sf = R(1).regressor_events{end};
% event_weights_sf = eye(numel(event_list_sf));                       
                        
chanlocs = R(1).EEG.chanlocs;
tResponse = R(1).tResponse{end};


%% Set up (get single-subj stats)
doPlot = false;
[contrastFns_sf, contrastVar_sf, contrastZ_sf]  = SetUpTopLevelGlm_flex(R(:,iSuffix),event_list_sf,event_weights_sf,[],[],doPlot);
%% Run (get group-level stats)
multcorrect = 'none';
[group_RF_sf,group_P_sf] = RunTopLevelGlm_EEG(contrastFns_sf,contrastVar_sf,multcorrect);
group_Z_sf = norminv(group_P_sf);

%% Plot group results as ERPs
% PlotResponseFnsForThesis({R(:,4)},{[1:4 8:11]},chansToPlot,iLevel_vec,addSquareRf);
legendstr_sf = {'T0-D0','D1-D0','T1-D0'};
% legendstr_sf = {'T0-D0','D1-D0','T1-D0','D2-D0','T2-D0'};
% legendstr_sf = event_list_sf;
chansToPlot = {'FZ';'CZ';'PZ';'OZ'};
Cmap = GetSquaresEventColormap(event_list_sf(2:4));
% Cmap = GetSquaresEventColormap(event_list_sf(2:6));

figure(61); clf;
PlotResponseFnsGrid(group_Z_sf, legendstr_sf,tResponse,chanlocs,chansToPlot,Cmap);
MakeFigureTitle(suff);
set(gcf,'Position',[1 268 568 1238]);

%% Plot scalp maps
tBinCenters = 37.5:75:750;
tBinWidth = 75;

figure(2); clf;
[sm_all,sm] = GetScalpMaps(group_RF_sf,tResponse,tBinCenters,tBinWidth);
PlotScalpMaps(sm_all,chanlocs,[],tBinCenters-tBinWidth/2,legendstr_sf);
% Annotate figure
MakeFigureTitle(sprintf('Level 3 results for %d sf Subjects',numel(R(:,iSuffix))));


%% Plot statistics
cthresh = 1.96; % z score for 2-tailed p=0.05

figure(3); clf;
group_Z_sf = norminv(group_P_sf);
[sm_all,sm] = GetScalpMaps(group_Z_sf,tResponse,tBinCenters,tBinWidth,cthresh);
PlotScalpMaps(sm_all,chanlocs,[],tBinCenters-tBinWidth/2,legendstr_sf);
% Annotate figure
MakeFigureTitle(sprintf('Level 3 results for %d sf Subjects',numel(R(:,iSuffix))));
