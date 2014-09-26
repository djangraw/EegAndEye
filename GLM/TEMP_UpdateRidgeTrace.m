for iSubj=1:numel(subjects)
    fprintf('%s Subject %d...\n',datestr(now,16),iSubj);
    cd(basedir)
    cd(folders{iSubj});
    fprintf('%s Loading...\n',datestr(now,16))
    R = load(sprintf('%s-%d-GLMresults-Type-v3pt6-RampUp',experiment,subjects(iSubj)));
    fprintf('%s Getting Matrices...\n',datestr(now,16))
    [~,Y,Xmean,Xrss] = GetGlmMatrices(R.EEG,R.regressor_events{1}(iEvents),R.influence{1},R.artifact_events,R.artifact_influence{1},demeanX,normalizeX);
    %%
    Ymean = zeros(1,size(Y,2)); % don't de-mean Y
    Yrss = nan(1,size(Y,2));
    for i=1:size(Y,2)
        foo = Y(:,i)-Ymean(i);
        Yrss(i) = sqrt(foo'*foo);
        Y(:,i) = foo/Yrss(i);
    end
    
    UpdateFile(sprintf('%s-%d-Type-v3pt6-RampUp-RidgeTrace',experiment,subjects(iSubj)), {'Xrss','Xmean','Yrss','Ymean'},{Xrss,Xmean,Yrss,Ymean});    
    fprintf('%s Done!\n',datestr(now,16));
end