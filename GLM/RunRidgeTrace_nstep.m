function [betas, iBest, sigmasq_all, sigmasq_elec] = RunRidgeTrace_nstep(X_cell,Y_full,lambdas,nFolds,useRandomFolds)

% Created 2/26/15 by DJ.

[n,D] = size(Y_full);
L = numel(lambdas);

% Get folds
[iTrain,iTest] = GetFolds(n,nFolds,useRandomFolds);
[betas, sigmasq_all, sigmasq_elec] = deal(cell(1,nFolds));
% for each fold        
for iFold = 1:nFolds
    fprintf('%s ---Running fold %d/%d...\n',datestr(now,16),iFold,nFolds);
    % crop down to training & testing sets
%     fprintf('%s Setting up...\n',datestr(now,16))
    [X_train, X_test] = deal(cell(1,numel(X_cell)));
    for iLevel = 1:numel(X_cell)
        X_train{iLevel} = X_cell{iLevel}(iTrain{iFold},:);
        X_test{iLevel} = X_cell{iLevel}(iTest{iFold},:);
    end
    Y_train = Y_full(iTrain{iFold},:);
    Y_test = Y_full(iTest{iFold},:);

    %%% CALL THIS FN RECURSIVELY
    % run inner loop
    % get folds
    % for each fold
    % find beta (ridge trace)
    % find error
    % end inner loop

    if numel(X_train)>1
        fprintf('---> Entering inner loop...\n')
        [betas_inner, iBest_inner, sigmasq_all_inner, sigmasq_elec_inner] = RunRidgeTrace_nstep(X_train(1:end-1),Y_train,lambdas,nFolds,useRandomFolds);
        % pick best xval betas                
        betas_best = mean(betas_inner{end}(:,:,iBest_inner(end),:),4); % mean across folds       
        % get residuals        
        Y_recon = X_test{end-1}*betas_best';
        Y_test = Y_test-Y_recon;
        Y_recon = X_train{end-1}*betas_best';
        Y_train = Y_train-Y_recon;
        fprintf('<--- Exiting inner loop...\n')
    else
        betas_inner = {};
        [iBest_inner,sigmasq_all_inner,sigmasq_elec_inner] = deal([]);
    end
    
    % find beta (ridge trace)
%     fprintf('%s Calculating betas...\n',datestr(now,16))
    betas_outer = RunRidgeTrace_simple(X_train{end},Y_train,lambdas); 
%     fprintf('%s Calculating Error...\n',datestr(now,16))
    % get MSE        
    sigmasq_all_outer = nan(L,1);
    sigmasq_elec_outer = nan(L,D);
    for j=1:L
        [sigmasq_all_outer(j), sigmasq_elec_outer(j,:)] = GetGlmMse(full(X_test{end}),Y_test,betas_outer(:,:,j));
%             fprintf('lambda = %g, o^2 = %g\n',lambdas(j), sigmasq_all(1,j,i));
    end
    [~, iBest_outer] = min(sigmasq_all_outer); % find best sigmasq value        
%     fprintf('%s Prepping results...\n',datestr(now,16))
    % combine inner/outer results
    betas = [betas_inner, {betas_outer}];
    sigmasq_all = cat(3,sigmasq_all_inner,sigmasq_all_outer);
    sigmasq_elec = cat(3,sigmasq_elec_inner,sigmasq_elec_outer);
    iBest = [iBest_inner, iBest_outer];    
end
fprintf('%s ---Done!\n',datestr(now,16))

% find beta (ridge trace)
% find error
% end outer loop

