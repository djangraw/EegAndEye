function [Zscore,Pval,Tstat,contrastFns,contrastVar] = GetGlmZscore(EEG0,results,multcompare,contrasts,doPlot)

% Get 'z scores' and residual power for a glm result.
%
% [Zscore,Pval,Tstat,contrastFns,contrastVar] = GetGlmZscore(EEG0,results,multcompare,contrasts,doPlot)
%
% INPUTS:
% - EEG0 is an eeglab dataset from a squares experiment, containing the 
% data from immediately BEFORE the GLM was run.  If left blank, results.EEG
% will be used.
% - results is the saved and reloaded output of RunEegGlm: a struct with
% fields 'regressor_events', 'responseFns', 'NewEEG', 'tResponse'.
% - multcompare is a string indicating the type of multiple-comparisons
% correction you want to use.  Supported: 'none' (default), 'fdr' (false
% discovery rate), 'bonferroni'.
% - contrasts is a PxM matrix, where P is the number of regressor events
% and M is the number of contrasts [default is eye(P)].
% - doPlot is a binary value indicating whether you'd like to plot the
% results in the current figure. [default is true]
%
% OUTPUTS:
% - Zscore is a matrix of the same size as results.responseFns, indicating
% the responseFn 'betas' divided by the residual power.  A larger-magnitude 
% z score indicates a more significant correlation with the data.
% - Pval is a matrix of the same size that has the corresponding p values
% for each point.
% - Tstat is a matrix of the same size that has the corresponding T
% statistic for each point.
% -contrastVar is a DxTxMxN matrix inndicating the variance of each 
% contrast function estimate.
% -contrastZ is the single-subject Z score for each contrast function 
% estimate. No multiple-comparisons correction is performed.
% 
% Created 1/13/12 by DJ.
% Updated 1/25/12 by DJ - plot z score instead of p val, added zero lines
% Updated 5/10/12 by DJ - comments
% Updated 6/27/12 by DJ - fixed StdErr bug
% Updated 7/24/12 by DJ - updated for RunGlmGui results
% Updated 7/25/12 by DJ - fixed StdErr bug AGAIN, changed StdErr scale
% Updated 8/9/12 by DJ - update to work with event_weights, rejectepoch
% Updated 3/22/13 by DJ - use asymmetric influence, GetGlmRegressors_v2p0
% Updated 4/8/13 by DJ - add contrasts input, correction factors at Bryan's
%  advice
% Updated 5/1/13 by DJ - responseFns in cells
% Updated 3/19/14 by DJ - tResponse, influence, artifact_influence in cells
% Updated 8/7/14 by DJ - added doPlot input.
% Updated 8/11/14 by DJ - added option to load X instead of calculating it.
% Updated 8/15/14 by DJ - added degrees of freedom calc for ridge regressn
% Updated 8/25/14 by DJ - added parfors to speed calculations

% Handle inputs
if nargin<3 || isempty(multcompare)
    multcompare = 'none';
end
if nargin<4 
    contrasts = [];
end
if nargin<5
    doPlot = true; % plot by default.
end

% Parse results struct
results = UpdateGlmResultsFormat(results);
% [regressor_events, responseFns, tResponse] = deal(0); % just to keep compiler warnings away
UnpackStruct(results); % create variables regressor_events, responseFns, EEG, tResponse
if ~isempty(EEG0) % unless the subject wanted to use the EEG data in the results struct
    EEG = EEG0;
end
lambda = results.lambda; % to satisfly parallel job
iLevel = results.iLevel; % to satisfly parallel job
clear EEG0 results

if ~exist('tResponse','var'), tResponse = []; end; % fool the compiler
if iscell(tResponse), tResponse = tResponse{iLevel}; end % extract time vector
% Set up
dt = mode(diff(tResponse)); % ms between samples
regressor_range = round(influence{iLevel}/dt); % how many samples should each response function extend?
artifact_range = round(artifact_influence{iLevel}/dt); % how many samples should each artifact affect?
Nh = size(responseFns{iLevel},2); % # time points in responseFn
dt = 1000/EEG.srate;
t = (1:EEG.pnts*EEG.trials)*dt; % time vector for EEG

% Find event times
Nr = size(responseFns{iLevel},3); % # regressors
D = EEG.nbchan;
p = Nr*Nh;

% Load X matrix if possible
iMiddle = strfind(filenames{iLevel},'-GLMresults-');
prefix = filenames{iLevel}(1:iMiddle-1);
suffix = filenames{iLevel}(iMiddle+length('-GLMresults-'):end);
if 0 %strcmp(suffix,'Type-v3pt5-RampUp.mat') && strncmp(prefix,'sf-',3)
    % This may or may not be an actual shortcut...
    tic
    disp('Loading X & Y matrices...');
    load(sprintf('%s-%s-GlmMats',prefix,suffix(1:end-4))); % loads X and Y
    
    % Get residuals
    r=Y; % residuals
    for i = 1:iLevel
        if ~isempty(responseFns{i})
            fprintf('Level %d...\n',i);
            foo = permute(responseFns{i},[2 3 1]);
            H = reshape(foo,[size(foo,1)*size(foo,2),size(foo,3)]);
            reconstructed = (X*H);
            
            r = reconstructed-r;
        end         
    end
    toc
    
else
    % Calculate the X matrix and find residuals
    tic
    disp('Recalculating X & Y matrices...');
    event_times = cell(1,Nr);
    event_weights = cell(1,Nr);
    for i=1:Nr
        if isfield(EEG.etc,'rejectepoch')
            isGoodEvent = strcmp(regressor_events{iLevel}{i},{EEG.event(:).type}) & ~EEG.etc.rejectepoch([EEG.event(:).epoch])';
        else
            isGoodEvent = strcmp(regressor_events{iLevel}{i},{EEG.event(:).type});
        end
        event_times{i} = [EEG.event(isGoodEvent).latency]*dt;
        if isfield(EEG.etc,'ureventweights') % v2pt1 results and after
            event_weights{i} = EEG.etc.ureventweights([EEG.event(isGoodEvent).urevent]);
        elseif isfield(EEG.etc,'eventweights') % v2pt0 results
            event_weights{i} = EEG.etc.eventweights(isGoodEvent);
        else % v1ptX results
            event_weights{i} = ones(size(event_times{i}));
        end
    end

    % Find blink events
    artifact_times = [EEG.event(ismember({EEG.event(:).type}, artifact_events)).latency]*dt;

    % Get regressors
    disp('Getting regressors...');
    % [~,S] = GetGlmRegressors(t,event_times,artifact_times,regressor_range,event_weights,stddev);
    [~,S] = GetGlmRegressors_v2p0(t,event_times,artifact_times,regressor_range,artifact_range,event_weights,stddev);

    % Get cropped regressor matrix for later calculations
    isInPlay = any(S~=0,2);
    X = full(S(isInPlay,:)); % cropped regressor matrix
    clear S
    % Get residuals after fit
%     disp('Getting residuals...');
%     for i = 1:iLevel
%         if ~isempty(responseFns{i})
%             fprintf('Level %d...\n',i);
%             EEG = SubtractOutGlmResponses(EEG,responseFns{i},influence{i},regressor_events{i},stddev);
%         end
%     end
%     % Get residuals
%     Eres = reshape(EEG.data,EEG.nbchan,EEG.pnts*EEG.trials);        
%     r = Eres(:,isInPlay)'; % residuals

    % For Ridge Regression: Center & standardize S matrix
    if strcmp(method,'ridge')                
        if demeanX
            fprintf('De-meaning and ');
            mnX = mean(X,1);
        else
            fprintf('NOT de-meaning, but ');
            mnX = zeros(1,size(X,2));
        end
        fprintf('Standardizing design matrix...\n')
        p = size(X,2);
        rssX = nan(size(mnX));
        for i=1:p  
            foo = (X(:,i) - mnX(i));
            rssX(i) = sqrt(sum(foo.^2)); % root sum square
            X(:,i) = foo / rssX(i);
        end
    end
        
    % Faster residual calculation
    disp('Getting residuals...');
    r = reshape(EEG.data,EEG.nbchan,EEG.pnts*EEG.trials)';        
    clear EEG;
    r = r(isInPlay,:);
    for i = 1:iLevel
        if ~isempty(responseFns{i})
            fprintf('Level %d...\n',i);            
            H = reshape(responseFns{i},[D p])';
            reconstructed = (X*H);
            
            r = r - reconstructed;
        end         
    end
    
    
    % clean up
    clear S Eres reconstructed isGoodEvent isInPlay event_times event_weights
    toc
end


% Get sum of squared residuals and covar matrix
sigmasq = sum(r.^2); % sum of squared residuals
XTX = X'*X; % covariance matrix
n = size(X,1); % number of samples

% Free up more memory!
clear EEG r X t;


% Calculate degrees of freedom and covariance matrix inverse.
if strcmp(method,'ridge')    
%     dof = n-trace(X * (X'*X + lambda*eye(p))^-1 * X'); % OVERLOADS MEMORY!!!

    E = eig(XTX);% get eigenvalues   
    
    if length(lambda{iLevel})>1
        X_inv = zeros(p,p,D);
        dof = zeros(1,D);
        parfor j=1:D        
            fprintf('%d/%d: Calculating dof... ',j,D);
            dof(j) = n-sum(E./(E+lambda{iLevel}(j)));  % df(lambda) = sum(i=1:p)(d_i^2/(d_i^2+lambda))
            fprintf('dof = %g. Calcluating inverse... ',dof(j));
            tic;
            X_inv(:,:,j) = (XTX + lambda{iLevel}(j)*eye(p))^(-1);
            t = toc;
            fprintf('Done (took %.1f seconds).\n',t);        
        end
    else
        fprintf('Calculating dof... ');
        dof = n-sum(E./(E+lambda{iLevel}));  % df(lambda) = sum(i=1:p)(d_i^2/(d_i^2+lambda))
        fprintf('dof = %g. Calcluating inverse... ',dof);
        tic;
        X_inv = (XTX + lambda{iLevel}*eye(p))^(-1);
        t = toc;
        dof = repmat(dof,1,D);
        fprintf('Done (took %.1f seconds).\n',t);        
    end
else
    fprintf('Calcluating inverse... ');
    tic;
    X_inv = XTX^(-1);
    t = toc;
    fprintf('Done (took %.1f seconds).\n',t);
    dof = n-p; % degrees of freedom
    fprintf('Measurement has %.1f degrees of freedom.\n',dof);
    dof = repmat(dof,1,D);
end

if isempty(contrasts)
    contrasts = eye(p);
end

% normalize contrasts
if strcmp(method,'ridge')  
    for i=1:p
        contrasts(i,:) = (contrasts(i,:)-mnX(i))/rssX(i);
    end
end

% correction_fact = contrasts'*X_inv*contrasts; % correction factor
% beta_var = (sigmasq/dof)*correction_fact;
% residStdErr = sqrt(beta_var);

Nc = size(contrasts,2)/Nh; % # contrasts
beta = reshape(responseFns{iLevel},[D p]);


fprintf('Getting correction factors...')
tic;
% Get correction factors
if size(X_inv,3)>1
    cfact_cell = cell(1,D);
    parfor j=1:D
        for i=1:(Nc*Nh)
            cfact_cell{j}(i) = contrasts(:,i)'*X_inv(:,:,j)*contrasts(:,i); % correction factor to account for linear dependence between contrasts
        end
    end
    correction_fact = cat(1,cfact_cell{:});
else
    correction_fact = nan(D,Nc*Nh);
    for i=1:(Nc*Nh)
        correction_fact(:,i) = contrasts(:,i)'*X_inv*contrasts(:,i); % correction factor to account for linear dependence between contrasts
    end
end
t=toc;
fprintf('Done! Took %.1f seconds.\n',t);

fprintf('Calculating contrasts...')
tic;
contrastFns = zeros(D,Nh,Nc); % like responseFns
contrastVar = zeros(D,Nh,Nc);
for j=1:D
    for i=1:(Nc*Nh)    
        % Get inidces
        iCF = ceil(i/Nh); % which contrast fn?
        iT = i-Nh*(iCF-1); % which time point in that CF?
           
%         contrastVar(:,iT,iCF) = sigmasq'*correction_fact(:,i)/dof; % variance of contrast function
        contrastVar(j,iT,iCF) = sigmasq(j)*correction_fact(j,i)/dof(j); % variance of contrast function
        
        % Get response fn of contrast    
%         contrastFns(:,iT,iCF) = beta*contrasts(:,i);
        contrastFns(j,iT,iCF) = beta(j,:)*contrasts(:,i); % convert response functions to contrast functions
    end    
end
t=toc;
fprintf('Done! Took %.1f seconds.\n',t);

% take sqrt to get std error (???) of residuals
residStdErr = sqrt(contrastVar); 

% Convert to T statistics and P values.
Tstat = contrastFns ./ residStdErr;


%--------------------------------%
% Pval = 0; Zscore = 0; return; % TAKE THIS OUT!!!!!!
%--------------------------------%
Pval_start = nan(size(Tstat));
for j=1:D
    Pval_start(j,:) = tcdf(Tstat(j,:),dof(j));
end
% Apply multiple-comparisons correction
switch multcompare
    case 'none'
        Pval = Pval_start;
    case 'fdr' % false discovery rate
        isHigh = Pval_start>0.5;
        Pval_start(isHigh) = 1-Pval_start(isHigh);        
        Pval = reshape(mafdr(Pval_start(:),'bhfdr',true),size(Pval_start));
        Pval(isHigh) = 1-Pval(isHigh);
    case 'bonferroni' % bonferroni correction
        isHigh = Pval_start>0.5;
        Pval_start(isHigh) = 1-Pval_start(isHigh);
        Pval = bonf_holm(Pval_start);
        Pval(Pval>0.5) = 0.5;
        Pval(isHigh) = 1-Pval(isHigh);
    otherwise
        error('multiple comparisons method not recognized!');
end

% Convert to Z scores
Zscore = norminv(Pval);

% Plot results
if doPlot
    disp('Plotting results...');
    clf;

    for i=1:Nc
        subplot(3,Nc,i); hold on;
        plot(tResponse, contrastFns(:,:,i)');
        plot(tResponse, mean(contrastFns(:,:,i),1),'k','linewidth',2);
        plot([0 0],[-6 6],'k');
        plot([tResponse(1),tResponse(end)],[0 0],'k');
    %     title(regressor_events{iLevel}{i})
        ylabel('Response fn')
        ylim([-6 6])
        subplot(3,Nc,Nc+i); hold on;
        plot(tResponse, residStdErr(:,:,i)');
        plot(tResponse, mean(residStdErr(:,:,i),1),'k','linewidth',2);
        ylabel('StdErr of residuals');
        ylim([0 1])
        plot([0 0],get(gca,'ylim'),'k');
        plot([tResponse(1),tResponse(end)],[0 0],'k');
        subplot(3,Nc,2*Nc+i); hold on;   
        plot(tResponse, Zscore(:,:,i)');
        plot(tResponse, mean(Zscore(:,:,i),1),'k','linewidth',2);
        ylabel('Z score')
        plot([tResponse(1),tResponse(end)],[0 0],'k');
        plot([tResponse(1),tResponse(end)],norminv([.025 .025]),'k--');
        plot([tResponse(1),tResponse(end)],norminv([.975 .975]),'k--');
        plot([0 0],[-4 4],'k');
        ylim([-4 4])

        xlabel('time from event (ms)')    


    end
end
disp('Success!');