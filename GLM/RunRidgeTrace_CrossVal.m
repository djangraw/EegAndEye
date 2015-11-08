function RunRidgeTrace_CrossVal(file_in,file_out,lambdas,nFolds,params)

% RunRidgeTrace_CrossVal(filename,lambdas,nFolds,params)
%
% INPUTS:
%
% OUTPUTS:
%
% Created 9/11/14 by DJ.
% Updated 2/3/15 by DJ - added useRandomFolds option, stopped results print
% Updated 2/24/15 by DJ - fixed GetGlmMse input bug!

UnpackStruct(params);
if ~exist('useRandomFolds','var')
    useRandomFolds = false;
end
if ~exist('nRandomizations','var')
    nRandomizations = 1;
end
% Set up
fprintf('%s Loading...\n',datestr(now,16))
load(file_in); % loads X,Y,Xmean,Ymean,Xrss,Yrss

[n,p] = size(X);
D = size(Y,2);

%% split trials into folds
iTest = cell(nFolds,nRandomizations);
iTrain = cell(nFolds,nRandomizations);
folds = ceil(linspace(eps,nFolds,n));
for j=1:nRandomizations
    if useRandomFolds    
        thisFolds = folds(randperm(numel(folds)));
    end
    for i=1:nFolds
        if useRandomFolds
            iTest{i,j} = find(thisFolds==i);
        else
            iTest{i,j} = round((i-1)*(n/nFolds))+1 : round(i * n/nFolds);
        end
        iTrain{i,j} = setdiff(1:n,iTest{i,j});
    end
end

%%
% set up
sigmasq_all = nan(1,numel(lambdas),nFolds,nRandomizations);
sigmasq_elec = nan(numel(lambdas),D,nFolds,nRandomizations);
betas = nan(D,p,numel(lambdas),nFolds,nRandomizations);

% calculate
for k=1:nRandomizations
    for i=1:nFolds  
        if useRandomFolds
            fprintf('%s - randomization %d/%d, fold %d/%d...\n',datestr(now,16),k,nRandomizations,i,nFolds);    
        else
            fprintf('%s - fold %d/%d...\n',datestr(now,16),i,nFolds);    
        end
        Xtrain = X(iTrain{i,k},:);
        Xtest = X(iTest{i,k},:);
        Ytrain = Y(iTrain{i,k},:);
        Ytest = Y(iTest{i,k},:);
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
            betas(:,:,j,i,k) = ( (XTX+lambdas(j)*eye(p))^(-1)*XTY )'; % size Dxp
            [sigmasq_all(1,j,i,k), sigmasq_elec(j,:,i,k)] = GetGlmMse(full(Xtest),Ytest,betas(:,:,j,i,k));
%             fprintf('lambda = %g, o^2 = %g\n',lambdas(j), sigmasq_all(1,j,i));
        end
    end
end

% save results
save(file_out,'sigmasq_all','sigmasq_elec','betas','iTest','iTrain','lambdas');

