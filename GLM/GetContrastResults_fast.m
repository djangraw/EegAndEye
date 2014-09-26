% GetContrastResults_fast
%
% Produce group stats based on results made by GetRidgeTrace.
%
% Created 8/28/14 by DJ.
%% Load results
% experiment = 'sq';
% subjects = [9:11, 13:15 17:27];
experiment = 'sf';
subjects = [1:10 12:13];
suffix = 'Type-v3pt6-RampUp-LOSO';
% experiment = 'sf3';
% subjects = 1:12;
% suffix = 'Type-v3pt6-RampUp-RidgeTrace';
clear R
for iSubj=1:numel(subjects)
    fprintf('Subject %d...\n',iSubj);
    R(iSubj) = load(sprintf('%s-%d-%s',experiment,subjects(iSubj),suffix));    
end

%% Extract variables
% 
% % DECLARE LAMBDA
% lambda = 0.2;
iLambda = find(R(1).lambdas==lambda);

[XTX_all,sigmasq,n,beta_all,Xmean,Xrss] = deal([]);
for i=1:numel(R)
    XTX_all = cat(3,XTX_all,R(i).XTX);
    sigmasq = cat(1,sigmasq,R(i).sigmasq_all(iLambda));
    n = cat(1,n,R(i).n);
    beta_all = cat(3,beta_all,R(i).betas(:,:,iLambda));
    Xmean = cat(1,Xmean,R(i).Xmean);
    Xrss = cat(1,Xrss,R(i).Xrss);
end

%% Get contrasts
R0 = load(sprintf('%s-%d-GLMresults-Type-v3pt6-RampUp',experiment,subjects(1)));
if strcmp(experiment,'sf3')
    iEvents = 2:25;
else
    iEvents = 2:19;
end
event_list = R0.regressor_events{end}(iEvents);
tResponse = R0.tResponse{end};
chanlocs = R0.EEG.chanlocs;
% rule = 'T0vD0';
[contrasts0,legendstr,sequence,tContrast] = GetSequenceContrast(experiment,rule,event_list,tResponse);

%% Get single-subject contrasts

% set up contrast functions
N = numel(R);
[p,D]= size(R(1).XTY);
T = length(tResponse);
M = size(contrasts0,2)/T; % number of contrasts
[contrastFns, contrastVar, contrastZ] = deal(nan(D,T,M,N));
multcorrect = 'none';
% find contrast functions
for iSubj=1:N
    fprintf('---subject %d/%d...\n',iSubj,N);
    XTX = XTX_all(:,:,iSubj);
    % Calculate and plot
%     [contrastZ(:,:,:,i),~,~,contrastFns(:,:,:,i),contrastVar(:,:,:,i)] = GetGlmZscore(R(i).EEG,R(i),multcorrect,contrasts,doPlot);
    fprintf('Calculating dof... ');    
    E = eig(XTX);
    dof = n(iSubj)-sum(E./(E+lambda));  % df(lambda) = sum(i=1:p)(d_i^2/(d_i^2+lambda))
    fprintf('dof = %g. \nCalcluating inverse... ',dof);
    tic;
    X_inv = (XTX + lambda*eye(p))^(-1);
    t = toc;
    
    
    contrasts = contrasts0;
    for i=1:p
        contrasts(i,:) = (contrasts(i,:)-Xmean(iSubj,i))/Xrss(iSubj,i);
    end


    beta = beta_all(:,:,iSubj);


    fprintf('Getting correction factors...')
    tic;
    % Get correction factors
    correction_fact = nan(D,M*T);
    for i=1:(M*T)
        correction_fact(:,i) = contrasts(:,i)'*X_inv*contrasts(:,i); % correction factor to account for linear dependence between contrasts
    end
    t=toc;
    fprintf('Done! Took %.1f seconds.\n',t);

    fprintf('Calculating contrasts...')
    tic;
%     contrastFns = zeros(D,T,M); % like responseFns
%     contrastVar = zeros(D,T,M);
    for j=1:D
        for i=1:size(contrasts,2)    
            % Get inidces
            iCF = ceil(i/T); % which contrast fn?
            iT = i-T*(iCF-1); % which time point in that CF?

    %         contrastVar(:,iT,iCF) = sigmasq'*correction_fact(:,i)/dof; % variance of contrast function
            contrastVar(j,iT,iCF,iSubj) = sigmasq(iSubj)*correction_fact(j,i)/dof; % variance of contrast function

            % Get response fn of contrast    
    %         contrastFns(:,iT,iCF) = beta*contrasts(:,i);
            contrastFns(j,iT,iCF,iSubj) = beta(j,:)*contrasts(:,i); % convert response functions to contrast functions
        end    
    end
    t=toc;
    fprintf('Done! Took %.1f seconds.\n',t);
    
    disp('Getting statistics...')
    tic
    % take sqrt to get std error (???) of residuals
    residStdErr = sqrt(contrastVar); 

    % Convert to T statistics and P values.
    Tstat = contrastFns ./ residStdErr;

    Pval = nan(size(Tstat));
    for j=1:D
        Pval(j,:) = tcdf(Tstat(j,:),dof);
    end
    contrastZ(:,:,:,iSubj) = norminv(Pval);
    
    t=toc;
    fprintf('Done! Took %.1f seconds.\n',t);
      
end

%% Plot single-subject results

if doPlot
    %% Plot responses
    subjstr = cell(1,N);
    for i=1:N
        subjstr{i} = sprintf('Subj %d',subjects(i));
    end
    figure(201); clf;
    [sm_all,sm] = GetScalpMaps(permute(contrastFns,[1 2 4 3]),tResponse,tBinCenters,tBinWidth);
    PlotScalpMaps(sm_all,chanlocs,[],tBinCenters-tBinWidth/2,subjstr);
    MakeFigureTitle(sprintf('%s %s, %s contrast, Single-Subject RFs',experiment,suffix,rule));

    %% Plot Z scores
    cthresh = 1.96; % z score for 2-tailed p=0.05
    figure(202); clf;
    [sm_all,sm] = GetScalpMaps(permute(contrastZ,[1 2 4 3]),tResponse,tBinCenters,tBinWidth,cthresh);
    PlotScalpMaps(sm_all,chanlocs,[],tBinCenters-tBinWidth/2,subjstr);
    MakeFigureTitle(sprintf('%s %s, %s contrast, Single-Subject Z scores',experiment,suffix,rule));
end

%% Get group-level stats
multcorrect = 'none';
[group_RF,group_P] = RunTopLevelGlm_EEG(contrastFns,contrastVar,multcorrect);
group_Z = norminv(group_P);


%% Plot group results


if doPlot
    %% Plot group results as ERPs
    chansToPlot = {'FZ';'CZ';'PZ';'OZ'};
    % Cmap = GetSquaresEventColormap({'pT_{0/2}'});
    % Cmap = GetSquaresEventColormap({'pT_{0/2}','pD_{1/2}','pT^*_{1/2}'});
    Cmap = GetSquaresEventColormap({'pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pT^*_{2/3}'});
    % Cmap = GetSquaresEventColormap(event_list_sf(2:6));
    Cmap = {'r','g','b','c','m','y','k'};
    figure(61); clf;
    PlotResponseFnsGrid(group_RF, legendstr,tResponse,chanlocs,chansToPlot,Cmap);
    % MakeFigureTitle(suff);
    set(gcf,'Position',[1 268 568 1238]);
    for i=1:4
    subplot(4,1,i);
    title('')
    ylabel(chansToPlot{i});
    end
    subplot(4,1,1);
    title(sprintf('%s %s, Group RFs',experiment,suffix))

    figure(71); clf;
    PlotResponseFnsGrid(group_Z, legendstr,tResponse,chanlocs,chansToPlot,Cmap);
    % MakeFigureTitle(suff);
    set(gcf,'Position',[569 268 568 1238]);
    for i=1:4
    subplot(4,1,i);
    title('')
    ylabel(chansToPlot{i});
    end
    subplot(4,1,1);
    title(sprintf('%s %s, Group Z scores',experiment,suffix))



    %% Plot scalp maps
    tBinCenters = 37.5:75:1000;%750;
    tBinWidth = 75;

    figure(2); clf;
    [sm_all,sm] = GetScalpMaps(group_RF,tResponse,tBinCenters,tBinWidth);
    PlotScalpMaps(sm_all,chanlocs,[],tBinCenters-tBinWidth/2,legendstr);
    % Annotate figure
    MakeFigureTitle(sprintf('Group RFs for %d %s Subjects: %s',size(contrastFns,4),experiment,legendstr{1}));
    % set(gcf,'Position',[825 1352 1261 145]);


    % Plot statistics
    cthresh = 1.96; % z score for 2-tailed p=0.05
    
    figure(3); clf;
    % group_Z = norminv(group_P);
    [sm_all,sm] = GetScalpMaps(group_Z,tResponse,tBinCenters,tBinWidth,cthresh);
    PlotScalpMaps(sm_all,chanlocs,[],tBinCenters-tBinWidth/2,legendstr);
    % Annotate figure
    MakeFigureTitle(sprintf('Group Z scores for %d %s Subjects: %s',size(contrastFns,4),experiment,legendstr{1}));
    % set(gcf,'Position',[825 1202 1261 145]);
end

%% Save results

suffix = sprintf('-Type-v3pt6-RampUp-ridgep%02d',lambda*100);
% suffix = sprintf('-Type-v3pt6-RampUp-LosoLambda',lambda*100);
save(sprintf('%s%sResults-%scontrast',experiment,suffix,rule),'event_*','contrast*',...
    'group_*','chanlocs','tResponse','experiment','subjects','suffix',...
    'rule','sequence','tContrast','legendstr')


%% Reconstruct multiple contrasts

[group_RF,group_Z,group_P,contrastFns] = deal([]);
legendstr = {};
for i=1:numel(lambdas)
    lambda = lambdas(i);
    suffix = sprintf('-Type-v3pt6-RampUp-ridgep%02d',round(lambda*100));

    foo = load(sprintf('%s%sResults-%scontrast',experiment,suffix,rule));
    group_RF = cat(3,group_RF,foo.group_RF);
    group_Z = cat(3,group_Z,foo.group_Z);
    group_P = cat(3,group_P,foo.group_P);
    contrastFns = cat(3,contrastFns,foo.contrastFns);
%     legendstr = cat(2,legendstr,foo.legendstr);
    legendstr = cat(2,legendstr,sprintf('k=%.2f',lambda));
end

%% Correct group-level stats
multcorrect = 'fdr';
clear Pval_end
for i=1:size(group_P,3);
    Pval_start = group_P(:,:,i); 
    if all(isnan(Pval_start(:)))
        Pval_end(:,:,i) = Pval_start;
    else
        switch multcorrect
            case 'fdr' % false discovery rate
                isHigh = Pval_start>0.5;
                Pval_start(isHigh) = 1-Pval_start(isHigh);        
                Pval = reshape(mafdr(Pval_start(:),'bhfdr',true),size(Pval_start));
                Pval(isHigh) = 1-Pval(isHigh);
        end
        Pval_end(:,:,i) = Pval;
    end
end

group_Z = norminv(Pval_end);
