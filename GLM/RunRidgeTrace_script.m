% Created 9/11/14 by DJ.

% Use SaveOutGlmMatrices_script to get expts/subjects/folders

% Set parameters
lambdas = 0:0.04:1;
nFolds = 10;
params = struct('demeanX',0,'demeanY',0,'normailzeX',1,'normalizeY',1);

experiments = {'sq','sf','sf3'};

% suffix_in = 'Type-v3pt6-RampUp-Matrices';
% suffix_out = 'Type-v3pt6-RampUp';

% suffix_in = 'Type-v3pt6-Matrices';
% suffix_out = 'Type-v3pt6';

% suffix_in = 'Type-v3pt6-Peak-Matrices';
% suffix_out = 'Type-v3pt6-Peak';

%%
for i=1:numel(experiments)
    experiment = experiments{i};
    [subjects,basedir,folders] = GetSquaresSubjects(experiment);    
    for iSubj=1:numel(subjects)
        fprintf('%s %s Subject %d/%d...\n',datestr(now,16),experiment,iSubj,numel(subjects));
        cd(basedir)
        cd(folders{iSubj});
        file_in = sprintf('%s-%d-%s',experiment,subjects(iSubj),suffix_in);
        file_out = sprintf('%s-%d-%s-%dfold',experiment,subjects(iSubj),suffix_out,nFolds);
        if exist([file_out '.mat'],'file')
            disp('Skipping!');
        else
            RunRidgeTrace_CrossVal(file_in,file_out,lambdas,nFolds,params);
        end
    end
end

%% Get results
lambda_best = zeros(1,numel(experiments));
sigmasq_best = zeros(1,numel(experiments));
sigmasq_mean = zeros(numel(experiments),numel(lambdas));
sigmasq_foldmean = cell(1,numel(experiments));

for i=1:numel(experiments)
    experiment = experiments{i};
    subjects = subjects_cell{i};
    N = numel(subjects);
    sigmasq_foldmean{i} = zeros(N,numel(lambdas));
    for iSubj=1:N
        fprintf('%d/%d...\n',iSubj,N);
        file_out = sprintf('%s-%d-%s-%dfold',experiment,subjects(iSubj),suffix_out,nFolds);
        R = load(file_out);
        sigmasq_foldmean{i}(iSubj,:) = mean(R.sigmasq_all,3);
    end
    sigmasq_mean(i,:) = mean(sigmasq_foldmean{i},1);
    [~,iMin] = min(sigmasq_mean(i,:));
    lambda_best(i) = lambdas(iMin);
    sigmasq_best(i) = sigmasq_mean(i,iMin);
    fprintf('%s: %.2f\n',experiment,lambda_best(i));    
end 


%% Plot results
colors = 'rgb';
figure;
for i=1:numel(experiments)
    % set up plot
    subplot(1,numel(experiments),i);
    cla; hold on;
    % plot
    plot(lambdas,sigmasq_mean(i,:),colors(i));
    plot(lambda_best(i),sigmasq_best(i),[colors(i) '*']);    
    % annotate
    ylim([min(sigmasq_mean(i,:)),max(sigmasq_mean(i,2:end))])
    title(experiments{i})
    xlabel('bias param')
    ylabel('10-fold MSE')
end

%% Use Selected Lambda Values
[betas,sigmasq_all,sigmasq_elec] = deal(cell(1,numel(experiments)));
% D = 69;

for i=1:numel(experiments)
    experiment = experiments{i};
    subjects = GetSquaresSubjects(experiment);
    file_in = sprintf('%s-%d-%s',experiment,subjects(1),suffix_in);
    R = load(file_in);
    p = size(R.X,2);
    D = size(R.Y,2);
%     if strcmp(experiment,'sf3')
%         p = 2424;
%     else
%         p = 1818;
%     end
    N = numel(subjects);
    betas{i} = zeros(D,p,N);    
    sigmasq_all{i} = zeros(1,N);
    sigmasq_elec{i} = zeros(N,D);
    for iSubj=1:N
        fprintf('%d/%d...\n',iSubj,N);    
        file_in = sprintf('%s-%d-%s',experiment,subjects(iSubj),suffix_in);
        R = load(file_in);
        XTX = full(R.X'*R.X);
        XTY = full(R.X')*R.Y;
        betas{i}(:,:,iSubj) = ((XTX + lambda_best(i)*eye(p))^(-1)*XTY)';        
        
        [sigmasq_all{i}(iSubj), sigmasq_elec{i}(iSubj,:)] = GetGlmMse(R.X,R.Y,betas{i}(:,:,iSubj));
    end

end



%% Get group results
% old_suffix = 'GLMresults-Type-v3pt6-RampUp';
% suffix_in = 'Type-v3pt6-RampUp-Matrices';
% suffix_out = 'Type-v3pt6-RampUp';
old_suffix = 'GLMresults-Type-v3pt6-Peak';
suffix_in = 'Type-v3pt6-Matrices';
suffix_out = 'Type-v3pt6';

% rules = {'T0vD0','T1vT0','D1vD0','T1vD0','T+vD0','D+vD0','T*vD*','D*vD0'};
rules = {'D2vD1'};
multcorrect = 'none';
for iRule = 1:numel(rules)
    rule = rules{iRule};
    for iExp=3%1:numel(experiments)
        experiment = experiments{iExp};    
        [subjects,basedir,folders] = GetSquaresSubjects(experiment);
        cd(basedir);
        fprintf('--- %s %s: Setting up...\n',datestr(now,16),experiment);
        % set up contrast
        R = load(sprintf('%s-%d-%s',experiment,subjects(1),old_suffix));
        if strcmp(experiment,'sf3')
%             iEvents = 2:25; % for RampUp or Peak
            iEvents = 2:14; % for regular v3pt6
        else
%             iEvents = 2:19; % for RampUp or Peak
            iEvents = 2:12; % for regular v3pt6
        end
        p = length(R.tResponse{end})*numel(iEvents);
        event_list = R.regressor_events{end}(iEvents);
        tResponse = R.tResponse{end};
        chanlocs = R.EEG.chanlocs;
        % rule = 'T0vD0';
        [contrasts0,legendstr,sequence,tContrast] = GetSequenceContrast(experiment,rule,event_list,tResponse);


        %--- Get Single-Subject Contrasts
        lambda = lambda_best(iExp);
        % set up contrast functions
        D = numel(chanlocs);
        N = numel(subjects);    
        T = length(tResponse);
        M = size(contrasts0,2)/T; % number of contrasts
        [contrastFns, contrastVar, contrastZ] = deal(nan(D,T,M,N));
        multcorrect = 'none';
        % find contrast functions
        for iSubj=1:N
            fprintf('---subject %d/%d...\n',iSubj,N); 
            file_in = sprintf('%s-%d-%s',experiment,subjects(iSubj),suffix_in);
            R = load(file_in);
            XTX = full(R.X'*R.X);
            XTY = full(R.X')*R.Y;
            n = size(R.X,1);
            Xmean = R.Xmean;
            Xrss = R.Xrss;
            clear R;

            % Calculate and plot
        %     [contrastZ(:,:,:,i),~,~,contrastFns(:,:,:,i),contrastVar(:,:,:,i)] = GetGlmZscore(R(i).EEG,R(i),multcorrect,contrasts,doPlot);
            fprintf('Calculating dof... ');    
            E = eig(XTX);
            dof = n - sum(E./(E+lambda));  % df(lambda) = sum(i=1:p)(d_i^2/(d_i^2+lambda))
            fprintf('dof = %g. \nCalcluating inverse...\n',dof);
            tic;
            X_inv = (XTX + lambda*eye(p))^(-1);


            contrasts = contrasts0;
            for i=1:p
                contrasts(i,:) = (contrasts(i,:)-Xmean(i))/Xrss(i);
            end


            beta = betas{iExp}(:,:,iSubj);


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
                    contrastVar(j,iT,iCF,iSubj) = sigmasq_all{iExp}(iSubj)*correction_fact(j,i)/dof; % variance of contrast function

                    % Get response fn of contrast    
            %         contrastFns(:,iT,iCF) = beta*contrasts(:,i);
                    contrastFns(j,iT,iCF,iSubj) = beta(j,:)*contrasts(:,i); % convert response functions to contrast functions
                end    
            end
            t=toc;
            fprintf('Done! Took %.1f seconds.\n',t);
        end
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
        contrastZ = norminv(Pval);

        t=toc;
        fprintf('Done! Took %.1f seconds.\n',t);


        %% Get group-level stats
        multcorrect = 'none';
        [group_RF,group_P] = RunTopLevelGlm_EEG(contrastFns,contrastVar,multcorrect);
        group_Z = norminv(group_P);
        %%
        cd(basedir)
        save(sprintf('%s-%s-%dfoldResults-%scontrast.mat',experiment,suffix_out,nFolds,rule),'event_*','contrast*',...
            'group_*','chanlocs','tResponse','experiment','subjects','suffix_in','suffix_out',...
            'rule','sequence','tContrast','legendstr','lambda','multcorrect')



    end

end