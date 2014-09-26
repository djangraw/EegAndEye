% GetLosoLambdas_script.m
%
% Created 9/2/14 by DJ.

% experiment = 'sf3';
% subjects = 1:12;

% experiment = 'sf';
% subjects = [1:10 12:13];

experiment = 'sq';
subjects = [9:11, 13:15 17:27];


suffix = 'Type-v3pt6-RampUp-RidgeTrace';

lambda_all = nan(1,numel(subjects));
lambda_elec = nan(numel(subjects),D);

for iSubj = 1:numel(subjects)
    fprintf('%s Subject %d: Loading...\n',datestr(now,16),subjects(iSubj))
    R = load(sprintf('%s-%d-%s.mat',experiment,subjects(iSubj),suffix),'betas','sigmasq_all','lambdas','sigmasq_elec');
    
    fprintf('%s Getting Heuristic Estimates...\n',datestr(now,16))
    
    i0 = find(R.lambdas==0);
    beta0 = R.betas(:,:,i0);
    [D,p] = size(beta0);
    sigmasq_all0 = R.sigmasq_all(i0);
    lambda_all(iSubj) = (p*D)*sigmasq_all0/sum(beta0(:).^2);
    
    sigmasq_elec0 = R.sigmasq_elec(i0,:);
    for i=1:D
        lambda_elec(iSubj,i) = p*sigmasq_elec0(i)/sum(beta0(i,:).^2);
    end    
    
end

%% Calculate & Save LOSO lambdas
lambda_loso = nan(1,numel(subjects));
if any(isnan(lambda_all)) 
    error('Some lambdas are NaN!');
end
for iSubj = 1:numel(subjects)
    lambda_loso(iSubj) = mean(lambda_all([1:iSubj-1, iSubj+1:end]));
end
    
save(sprintf('%s-%s-lambdas-LOSO',experiment,suffix),'lambda_all','lambda_loso')




%% Perform Ridge Reg with loso lambdas

mat_suffix = 'Type-v3pt6-RampUp-Matrices';
out_suffix = 'Type-v3pt6-RampUp-LOSO';
cd(basedir)
suffix = 'Type-v3pt6-RampUp-RidgeTrace';
load(sprintf('%s-%s-lambdas-LOSO',experiment,suffix));
for iSubj = 1:numel(subjects)
    
    cd(basedir);
    cd(folders{iSubj});
    
    fprintf('--- %s Subject %d: Loading...\n',datestr(now,16),subjects(iSubj))
    load(sprintf('%s-%d-%s.mat',experiment,subjects(iSubj),mat_suffix));
    
    fprintf('--- %s Multiplying...\n',datestr(now,16))
    XTX = X'*X;
    XTY = X'*Y;
    n = size(X,1);
    
    fprintf('--- %s Solving using LOSO lambda...\n',datestr(now,16))
    lambda = lambda_loso(iSubj);
    betas = ( (XTX + lambda*eye(size(XTX,1)) )^-1 * XTY)';    
    [sigmasq_all, sigmasq_elec] = GetGlmMse(X,Y,betas);
    fprintf('lambda = %g, o^2 = %g\n',lambda, sigmasq_all);
    %% save    
    new_filename = sprintf('%s-%d-%s',experiment,subjects(iSubj),out_suffix);
    fprintf('%s Saving %s...\n',datestr(now,16),new_filename)
    save(new_filename,'XTX','XTY','lambda','sigmasq_all','sigmasq_elec','betas','n','Xmean','Xrss','Ymean','Yrss');
    fprintf('Done!\n')
end
cd(basedir)