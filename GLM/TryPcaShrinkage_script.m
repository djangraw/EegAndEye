
suffix = 'Type-v3pt6-RampUp-PCAshrink';
shrinkages = 0:200:1200;
%%
for iSubj=1:numel(subjects)
    fprintf('%s Subject %d...\n',datestr(now,16),iSubj);
    cd(basedir)
    cd(folders{iSubj});
    fprintf('%s Loading...\n',datestr(now,16))
    R = load(sprintf('%s-%d-GLMresults-Type-v3pt6-RampUp',experiment,subjects(iSubj)));
    fprintf('%s Getting Matrices...\n',datestr(now,16))
    [X,Y,Xmean,Xrss] = GetGlmMatrices(R.EEG,R.regressor_events{1}(iEvents),R.influence{1},R.artifact_events,R.artifact_influence{1},demeanX,normalizeX);
    
    %%
    Ymean = zeros(1,size(Y,2)); % don't de-mean Y
    Yrss = nan(1,size(Y,2));
    for i=1:size(Y,2)
        foo = Y(:,i)-Ymean(i);
        Yrss(i) = sqrt(foo'*foo);
        Y(:,i) = foo/Yrss(i);    
    end
    % clear memory
    clear R Xmean Xrss foo Ymean Yrss
    X = sparse(X); % to save memory
    
    %% Use PCs
    fprintf('%s Using PCA...\n',datestr(now,16))
    [betas] = cell(1,numel(shrinkages));
    p = size(X,2);
    
    [U,S,V] = svd(full(X));
    clear X Ymean Yrss
    
    fprintf('%s Solving...\n',datestr(now,16))
    for i=1:numel(shrinkages)
        fprintf('shrinkage = %d...\n',shrinkages(i));
%         [U,S,V] = svds(X,p-shrinkages(i));
%         X0 = U*S;
        X0 = U(:,1:p-shrinkages(i)) * S(1:p-shrinkages(i),1:p-shrinkages(i));
        
        XTX0 = X0'*X0;
        XTY0 = X0'*Y;
        beta0 = XTX0^-1*XTY0;
%         betas{i} = V*beta0;
        betas{i} = V(:,1:p-shrinkages(i))*beta0;
        
        
    end
    fprintf('%s Saving...\n',datestr(now,16));
    save(sprintf('%s-%d-%s.mat',experiment,subjects(iSubj),suffix),'shrinkages','betas','U','S','V');
    
    fprintf('%s Done!\n',datestr(now,16));
end