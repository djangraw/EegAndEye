% CompareGlmModelStats_script.m
%
% Compare 'partial' models to 'full' models with additional parameters.
% Use an F statistic to evaluate whether or not the additional parameters
% are worthwhile to include.
%
% Created 9/18/14 by DJ.

%% Get Stats
experiments = {'sq','sf','sf3'};
lambdas = 0:.04:1;

suffix_full = 'Type-v3pt6-Peak-10fold';
suffix_part = 'Type-v3pt6-10fold';
suffix_full_mat = 'Type-v3pt6-Peak-Matrices';
suffix_part_mat = 'Type-v3pt6-Matrices';

[lambda_best, sigmasq_best] = deal(zeros(2,3));
fprintf('=== Getting 10-fold lambdas ===\n')
for i=1:numel(experiments)    
    experiment = experiments{i};
    [subjects,basedir,folders] = GetSquaresSubjects(experiment);
    fprintf('--- %s: Starting %s Subjects ---\n',datestr(now,16),experiment)
    
    sigmasq_foldmean_full = zeros(numel(subjects),numel(lambdas));
    sigmasq_foldmean_part = zeros(numel(subjects),numel(lambdas));

    for iSubj = 1:numel(subjects)
        fprintf('%s: subject %d/%d...\n',datestr(now,16),iSubj,numel(subjects))
        subject = subjects(iSubj);
        cd(basedir);
        cd(folders{iSubj});
        
        % get sigmasq results - FULL
        R = load(sprintf('%s-%d-%s',experiment,subject,suffix_full));
        sigmasq_foldmean_full(iSubj,:) = mean(R.sigmasq_all,3);
        % get sigmasq results - PART
        R = load(sprintf('%s-%d-%s',experiment,subject,suffix_part));
        sigmasq_foldmean_part(iSubj,:) = mean(R.sigmasq_all,3);
        
    end
    
    sigmasq_mean = mean(sigmasq_foldmean_full,1);    
    [~,iMin] = min(sigmasq_mean);
    lambda_best(1,i) = lambdas(iMin);
    sigmasq_best(1,i) = sigmasq_mean(iMin);
    
    sigmasq_mean = mean(sigmasq_foldmean_part,1);    
    [~,iMin] = min(sigmasq_mean);
    lambda_best(2,i) = lambdas(iMin);
    sigmasq_best(2,i) = sigmasq_mean(iMin);
    fprintf('%s full: %.2f, part: %.2f\n',experiment,lambda_best(1,i),lambda_best(2,i));
end

%%

[betas,sigmasq_all,sigmasq_elec] = deal(cell(2,numel(experiments)));
fprintf('=== Getting betas and sigmasq values ===\n')
for i=1:numel(experiments)    
    experiment = experiments{i};
    [subjects,basedir,folders] = GetSquaresSubjects(experiment);
    fprintf('--- %s: Starting %s Subjects ---\n',datestr(now,16),experiment)

    R = load(sprintf('%s-%d-%s',experiment,subjects(1),suffix_full_mat));
    p_full = size(R.X,2);
    D = size(R.Y,2);
    R = load(sprintf('%s-%d-%s',experiment,subjects(1),suffix_part_mat));
    p_part = size(R.X,2);
    
    N = numel(subjects);
    betas{1,i} = zeros(D,p_full,N);    
    betas{2,i} = zeros(D,p_part,N);    
    [sigmasq_all{:,i}] = deal(zeros(1,N));
    [sigmasq_elec{:,i}] = deal(zeros(N,D));
    for iSubj = 1:numel(subjects)
        fprintf('%s: subject %d/%d...\n',datestr(now,16),iSubj,numel(subjects))
        subject = subjects(iSubj);

        % get matrices - FULL
        R = load(sprintf('%s-%d-%s',experiment,subject,suffix_full_mat));
        XTX = full(R.X'*R.X);
        XTY = full(R.X')*R.Y;
        % get betas
        betas{1,i}(:,:,iSubj) = ((XTX + lambda_best(1,i)*eye(size(XTX)))^(-1)*XTY)';        
        [sigmasq_all{1,i}(iSubj), sigmasq_elec{i}(iSubj,:)] = GetGlmMse(R.X,R.Y,betas{1,i}(:,:,iSubj));
        
        % Get matrices - PART
        R = load(sprintf('%s-%d-%s',experiment,subject,suffix_part_mat));
        XTX = full(R.X'*R.X);
        XTY = full(R.X')*R.Y;
        % get betas
        betas{2,i}(:,:,iSubj) = ((XTX + lambda_best(2,i)*eye(size(XTX)))^(-1)*XTY)';        
        [sigmasq_all{2,i}(iSubj), sigmasq_elec{i}(iSubj,:)] = GetGlmMse(R.X,R.Y,betas{2,i}(:,:,iSubj));
        
    end    
    
end

  

%% Get stats
[Fstat,Pval,sizes] = deal(cell(1,numel(experiments)));

for i=1:numel(experiments)    
    experiment = experiments{i};
    [subjects,basedir,folders] = GetSquaresSubjects(experiment);
    fprintf('--- %s: Starting %s Subjects ---\n',datestr(now,16),experiment)
    [Fstat{i},Pval{i}] = deal(zeros(1,numel(subjects)));    
    sizes{i} = repmat(struct('n',[],'D',[],'p',[],'k',[]),1,numel(subjects));
    for iSubj = 1:numel(subjects)
        fprintf('%s: subject %d/%d...\n',datestr(now,16),iSubj,numel(subjects))
        subject = subjects(iSubj);

        % get matrices - FULL
        R = load(sprintf('%s-%d-%s',experiment,subject,suffix_full_mat));
        [n,p] = size(R.X);
        D = size(R.Y,2);
        R = load(sprintf('%s-%d-%s',experiment,subject,suffix_part_mat));
        k = size(R.X,2);
        sizes{i}(iSubj) = struct('n',n,'D',D,'p',p,'k',k);
        
        SSE_full = sigmasq_all{1,i}(iSubj)*(n-p)*D;
        SSE_part = sigmasq_all{2,i}(iSubj)*(n-k)*D;
        
        Fstat{i}(iSubj) = ((SSE_part - SSE_full)/((p-k)*D)) / (SSE_full / ((n-p)*D));
        Pval{i}(iSubj) = 1-fcdf(Fstat{i}(iSubj),(p-k)*D,(n-p)*D);
    end
end
        
%% save out results
save(sprintf('%s_vs_%s-Betas',suffix_full,suffix_part),'suffix_full','suffix_part','lambda_best','betas','sigmasq_all','sigmasq_elec','Fstat','Pval','sizes');

%% Plot results
figure(456); clf;
for i=1:numel(experiments)
    subplot(numel(experiments),1,i);
    bar(Pval{i});
%     foo = cat(1,sigmasq_all{:,i});
%     foo = cat(1,R1.sigmasq_all{1,i}, R2.sigmasq_all{:,i});    
%     foo = R1.sigmasq_all{1,i}-R2.sigmasq_all{2,i};
%     bar(foo');
%     legend('RampUp','Peak','Reduced')
%     legend('Error(RampUp)-Error(Peak)')
    ylabel(experiments{i})       
end
xlabel('subject');
subplot(numel(experiments),1,1);
title('sigmasq(RampUp) - sigmasq(v3pt6)');
        
        