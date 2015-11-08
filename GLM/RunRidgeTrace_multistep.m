function [sigmasq_all,sigmasq_elec] = RunRidgeTrace_multistep(X_full,Y_full,lambdas,nFolds,useRandomFolds)

D = size(Y_full,2);
L = numel(lambdas);
% Get Xval
if iscell(X_full)
    n = nan(1,numle(X_full));
    [iTrain,iTest] = deal(cell(1,numel(X_full)));
    n(1) = size(Y_full,1);
    [iTrain{1},iTest{1}] = GetFolds(n(1),nFolds,useRandomFolds);
    for i=2:numel(X_full)        
        n(i) = numel(iTrain{i-1});
        [iTrain_new,iTest_new] = GetFolds(n(i),nFolds,useRandomFolds);
        iTrain{i} = iTrain{i-1}(iTrain_new);
        iTest{i} = iTrain{i-1}(iTest_new);
    end
    iTrain = fliplr(iTrain);
    iTest = fliplr(iTest);
    
    % train
    betas = cell(1,numel(X_full));
    [sigmasq_all, sigmasq_elec] = deal(cell(1,numel(X_full)));
    Yrss = zeros(numel(X_full),D);
    for i=1:numel(X_full) % level of analysis
        sigmasq_all{i} = nan(L,n(i));
        sigmasq_elec{i} = nan(L,D,n(i));
        betas{i} = nan(L,D,nFolds);
        for j=1:nFolds % fold
            Y = Y_full(iTrain{i}{j},:);        
            X = X_full{i}(iTrain{i}{j},:);
            % get betas
            [betas{i}(:,:,:,j),Yrss(i,:)] = RunRidgeTrace_simple(X,Y,lambdas);
            % get MSE
            Y = Y_full(iTest{i}{j},:);        
            X = X_full{i}(iTest{i}{j},:);
            for k=1:size(betas,3);    
                [err_all,err_elec] = GetGlmMse(full(X),Y,betas{i}(:,:,k,j));
                sigmasq_all{i}(k,iTest{i}{j}) = err_all;
                sigmasq_elec{i}(k,:,iTest{i}{j}) = err_elec;
            end            
            % subtract out
            [~,iBest] = min(sigmasq_all{i});
            betas_best{i} = mean(betas{i}(:,:,iBest,:),4);
        end
    end    
end
