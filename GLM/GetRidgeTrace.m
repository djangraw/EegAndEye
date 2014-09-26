experiment = 'sf3';
if strcmp(experiment,'sf3')
    iEvents = 2:25;
else
    iEvents = 2:19;
end
demeanX = 0;
normalizeX = 1;
lambdas = linspace(0,1,51);
% lambdas = linspace(0,.02,51);


old_suffix = 'GLMresults-Type-v3pt6-Peak';
new_suffix = 'Type-v3pt6-Peak-RidgeTrace';
% new_suffix = 'Type-v3pt6-RampUp-RidgeTrace1e-2';

switch experiment
    case 'sq'
        basedir = '/Users/dave/Documents/Data/Squares';
        subjects = [9:11, 13:15, 17:27];
        folders = {'2012-03-29-Pilot-S9', '2012-04-09-Pilot-S10', '2012-04-18-Pilot-S11', ...
            '2012-05-04-Pilot-S13', '2012-05-29-Pilot-S14', '2012-06-12-Pilot-S15', ...
            '2012-06-15-Pilot-S17', '2012-06-21-Pilot-S18', '2012-06-22-Pilot-S19',...
            '2012-06-25-Pilot-S20', '2013-04-15-Pilot-S21', '2013-04-16-Pilot-S22', ...
            '2014-02-13-Pilot-S23', '2014-02-17-Pilot-S24', '2014-02-20-Pilot-S25', ...
            '2014-03-04-Pilot-S26', '2014-03-10-Pilot-S27'};
    case 'sf'
        basedir = '/Users/dave/Documents/Data/SquaresFix';
        subjects = [1:10 12:13];
        folders = {'2013-03-05-Pilot-S1','2013-03-11-Pilot-S2','2013-03-18-Pilot-S3',...
            '2013-03-19-Pilot-S4','2013-04-01-Pilot-S5','2013-04-02-Pilot-S6',...
            '2013-04-09-Pilot-S7','2013-04-11-Pilot-S8','2013-04-11-Pilot-S9',...
            '2013-04-17-Pilot-S10','2014-02-17-Pilot-S12','2014-03-06-Pilot-S13'};
    case 'sf3'
        basedir = '/Users/dave/Documents/Data/SquaresFix3';
        subjects = 1:12;
        folders = {'2013-05-15-Pilot-S1'    '2013-12-02-Pilot-S2'    '2013-12-04-Pilot-S3'    '2013-12-04-Pilot-S4',...
            '2013-12-05-Pilot-S5'    '2013-12-16-Pilot-S6'    '2013-12-19-Pilot-S7'    '2014-02-04-Pilot-S8',...
            '2014-02-06-Pilot-S9'    '2014-02-06-Pilot-S10'   '2014-02-12-Pilot-S11'   '2014-02-18-Pilot-S12'};
end
%%

for iSubj=1:numel(subjects)
    fprintf('---Subject %d---\n',subjects(iSubj));
    % Check if we've done this subject already.
    cd(basedir);
    cd(folders{iSubj});
    if exist(sprintf('%s-%d-%s.mat',experiment,subjects(iSubj),new_suffix),'file')
        fprintf('Already Done... Skipping!\n');
        continue;
    end
    % Otherwise, load the data
    fprintf('%s Loading...\n',datestr(now,16))
    R = load(sprintf('%s-%d-GLMresults-%s',experiment,subjects(iSubj),old_suffix));
    fprintf('%s Getting Matrices...\n',datestr(now,16))
    [X,Y,Xmean,Xrss] = GetGlmMatrices(R.EEG,R.regressor_events{1}(iEvents),R.influence{1},R.artifact_events,R.artifact_influence{1},demeanX,normalizeX);

    %% Set up
    
%     X = Rnew.X(:,iEventsToUse);
%     % Normalize X
%     for i=1:size(X,2)
%         foo = X(:,i);%-mean(X(:,i));
%         X(:,i) = foo/sqrt(foo'*foo);
%     end

%     Y = double(Rnew.Y);
    % Normalize Y
    fprintf('Normalizing Y...\n');
%     Ymean = mean(Y,1); 
    Ymean = zeros(1,size(Y,2)); % don't de-mean Y~
    Yrss = nan(1,size(Y,2));
    for i=1:size(Y,2)
        foo = Y(:,i)-Ymean(i);
        Yrss(i) = sqrt(foo'*foo);
        Y(:,i) = foo/Yrss(i);
    end
    %% Add column of 1's
    % X = [X ones(size(X,1),1)];
    %%
    fprintf('%s Multiplying Matrices...\n',datestr(now,16))
    XTX = X'*X;
    XTY = X'*Y;
    %% ITERATE!
    fprintf('%s Getting Ridge Trace...\n',datestr(now,16))
    D = size(Y,2);
    p = size(X,2);
    n = size(X,1);
%     lambdas = linspace(0,1,51);
    sigmasq_all = nan(1,numel(lambdas));
    sigmasq_elec = nan(numel(lambdas),D);
    betas = nan(D,p,numel(lambdas));


%     parfor j=1:numel(lambdas)
%         betas(:,:,j) = ( (XTX+lambdas(j)*eye(p))^(-1)*XTY )'; % size Dxp    
%     end
%     fprintf('%s Getting Error Estimates...\n',datestr(now,16))
    for j=1:numel(lambdas)
        betas(:,:,j) = ( (XTX+lambdas(j)*eye(p))^(-1)*XTY )'; % size Dxp
        [sigmasq_all(j), sigmasq_elec(j,:)] = GetGlmMse(X,Y,betas(:,:,j));
        fprintf('lambda = %g, o^2 = %g\n',lambdas(j), sigmasq_all(j));
    end

    %% Get heuristic lambda estimate
    fprintf('%s Getting Heuristic Estimates...\n',datestr(now,16))
    beta0 = ( XTX^(-1)*XTY )';
    [sigmasq_all0, sigmasq_elec0] = GetGlmMse(X,Y,beta0);
    k_elec = zeros(1,D);
    for i=1:D
        k_elec(i) = p*sigmasq_elec0(i)/sum(beta0(i,:).^2);
    end
    k_all = (p*D)*sigmasq_all0/sum(beta0(:).^2);

    if doPlot
    %% Plot
        fprintf('%s Plotting Results...\n',datestr(now,16))
        figure(iSubj); clf;
        subplot(1,4,1);
        plot(lambdas,sigmasq_all);
        xlabel('bias parameter')
        ylabel('mean squared error')
        subplot(1,4,2);
        plot(lambdas,squeeze(betas(:,50,:))')
        xlabel('bias parameter')
        ylabel('beta estimates')
        fracsignchange = squeeze(mean(mean(abs(diff(betas>0,1,3)),1),2));
        subplot(1,4,3)
        plot(lambdas(1:end-1),fracsignchange);
        xlabel('bias parameter')
        ylabel('fraction of betas changing sign')
        subplot(1,4,4)
        topoplot(k_elec,chanlocs);
        colorbar
        title(sprintf('k estimates (k_{all} = %.3g)',k_all));
        MakeFigureTitle(R.filenames{1});
        set(gcf,'Position',[2 1246 1373 260]);
    end
    %% save    
    new_filename = sprintf('%s-%d-%s',experiment,subjects(iSubj),new_suffix);
    fprintf('%s Saving %s...\n',datestr(now,16),new_filename)
    save(new_filename,'XTX','XTY','lambdas','sigmasq_all','sigmasq_elec','betas','n','Xmean','Xrss','Ymean','Yrss','k_all','k_elec');
    fprintf('Done!\n')
end