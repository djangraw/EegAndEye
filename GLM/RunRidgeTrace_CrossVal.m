function RunRidgeTrace_CrossVal(file_in,file_out,lambdas,nFolds,params)

% RunRidgeTrace_CrossVal(filename,lambdas,nFolds,params)
%
% INPUTS:
%
% OUTPUTS:
%
% Created 9/11/14 by DJ.

UnpackStruct(params);

% Set up
fprintf('%s Loading...\n',datestr(now,16))
load(file_in); % loads X,Y,Xmean,Ymean,Xrss,Yrss

[n,p] = size(X);
D = size(Y,2);

%% split trials into folds
iTest = cell(1,nFolds);
iTrain = cell(1,nFolds);
for i=1:nFolds
    iTest{i} = round((i-1)*(n/nFolds))+1 : round(i * n/nFolds);
    iTrain{i} = setdiff(1:n,iTest{i});
end

%%
% set up
sigmasq_all = nan(1,numel(lambdas),nFolds);
sigmasq_elec = nan(numel(lambdas),D,nFolds);
betas = nan(D,p,numel(lambdas),nFolds);

% calculate
for i=1:nFolds  
    fprintf('%s - fold %d/%d...\n',datestr(now,16),i,nFolds);    
    Xtrain = X(iTrain{i},:);
    Xtest = X(iTest{i},:);
    Ytrain = Y(iTrain{i},:);
    Ytest = Y(iTest{i},:);
    % re-normalize
%     for k=1:p
%         normfactor = sqrt(Xtrain(:,k)'*Xtrain(:,k));
%         Xtrain(:,k) = Xtrain(:,k) / normfactor;
%         Xtest(:,k) = Xtest(:,k) / normfactor;
%     end
%     for k=1:D
%         normfactor = sqrt(Ytrain(:,k)'*Ytrain(:,k));
%         Ytrain(:,k) = Ytrain(:,k) / normfactor;
%         Ytest(:,k) = Ytest(:,k) / normfactor;
%     end

    
    
    XTX = full(Xtrain'*Xtrain);    
    XTY = full(Xtrain)'*Ytrain;
    % ridge trace
    for j=1:numel(lambdas)
        betas(:,:,j,i) = ( (XTX+lambdas(j)*eye(p))^(-1)*XTY )'; % size Dxp
        [sigmasq_all(1,j,i), sigmasq_elec(j,:,i)] = GetGlmMse(full(Xtest),Ytest,betas(:,:,j,i));
        fprintf('lambda = %g, o^2 = %g\n',lambdas(j), sigmasq_all(1,j,i));
    end
end

% save results
save(file_out,'sigmasq_all','sigmasq_elec','betas','iTest','iTrain','lambdas');

