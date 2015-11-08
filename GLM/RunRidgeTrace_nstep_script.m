function RunRidgeTrace_nstep_script()

% Created 2/27/15 by DJ.

% Use SaveOutGlmMatrices_script to get expts/subjects/folders

% Set parameters
lambdas = 0:0.04:1;
nFolds = 20;
params = struct('demeanX',0,'demeanY',0,'normailzeX',1,'normalizeY',1);

% experiments = {'sq','sf','sf3'};
experiments = {'sf'};

% suffixes_in = {'Square-v3pt6-Matrices', 'Type-v3pt6-Matrices', 'SqNum-v3pt6-Matrices', 'TargDis-v3pt6-Matrices'};
suffixes_in_cell = {{'TargDis-v3pt6-Matrices', 'Type-v3pt6-Matrices'}};
events_in_cell = {{[1:4], [3:11]}};
% suffixes_out = {'Square-v3pt6-randfold', 'Type-v3pt6-randfold', 'SqNum-v3pt6-randfold', 'TargDis-v3pt6-randfold'};
% suffixes_out = {'TargDis-v3pt6-randfold', 'Type-v3pt6-randfold'};
suffixes_out = {'TargDis-Type-v3pt6-randfold'};

useRandomFolds = true;
% nRandomizations = 1;

%%
for i=1:numel(suffixes_out)
    suffix_out = suffixes_out{i};
    suffixes_in = suffixes_in_cell{i};        
    fprintf('%s === %s ===\n',datestr(now,16),suffix_out);
    for iExp=1:numel(experiments)
        experiment = experiments{iExp};
        [subjects,basedir,folders] = GetSquaresSubjects(experiment);    
        for iSubj=12%1:numel(subjects)
            fprintf('%s %s Subject %d/%d...\n',datestr(now,16),experiment,iSubj,numel(subjects));
            cd(basedir)
            cd(folders{iSubj});            
            file_out = sprintf('%s-%d-%s-%dfold',experiment,subjects(iSubj),suffix_out,nFolds);
            if exist([file_out '.mat'],'file')
                disp('Skipping!');
            else
                X = cell(size(suffixes_in));
                for k=1:numel(suffixes_in)
                    file_in = sprintf('%s-%d-%s',experiment,subjects(iSubj),suffixes_in{k});
                    Mats=load(file_in);
                    X{k} = Mats.X;
                end
                Y = Mats.Y;                
                [betas, iBest, sigmasq_all, sigmasq_elec] = RunRidgeTrace_nstep(X,Y,lambdas,nFolds,useRandomFolds);
                save(file_out,'betas','iBest','sigmasq_all','sigmasq_elec','lambdas','nFolds','useRandomFolds','suffixes_in')
            end
        end
    end
end
fprintf('%s === Done!\n',datestr(now,16));

%% Plot results
if ~exist('chanlocs','var')
    EEG = pop_loadset('sf-1-all-40Hz-fs100-interp-noeog.set');
    chanlocs = EEG.chanlocs;
    iCz = find(strcmpi('CZ',{chanlocs.labels}));
    clear EEG;
end
tRF = 0:10:1000;
tMaps = 50:100:1000;
for iSubj = 12%1:numel(subjects)
    fprintf('%s Plotting %s Subject %d/%d...\n',datestr(now,16),experiment,iSubj,numel(subjects));
    file_out = sprintf('%s-%d-%s-%dfold',experiment,subjects(iSubj),suffix_out,nFolds);
    load(file_out); % betas, iBest,sigmasq_all,sigmasq_elec,lambdas,nFolds,useRandomFolds,suffixes_in    
    for k=1%:numel(suffixes_in)   
        file_in = sprintf('%s-%d-%s',experiment,subjects(iSubj),suffixes_in{k});
        Mats = load(file_in,'events');
        M = numel(Mats.events);
        D = size(betas{k},1);
        T = size(betas{k},2)/M;
        RF = reshape(betas{k}(:,:,iBest(k)),D,T,M);
        sm = GetScalpMaps(RF,tRF,tMaps,100);
        figure(10*iSubj + k);
        PlotScalpMaps(sm,chanlocs,[],tMaps,Mats.events);
    end
    figure(10*iSubj + numel(suffixes_in)+1);    
    plot(lambdas,squeeze(sigmasq_all));
    hold on
    PlotHorizontalLines(mean(Y(:).^2),'k:');
    xlabel('Lambda');
    ylabel('MSE');
    legend([suffixes_in, {'No fit'}])
    % subplots of channel errors
    figure(10*iSubj + numel(suffixes_in)+2);    
    foo = permute(sigmasq_elec,[2 1 3]);
    sm = GetScalpMaps(foo,lambdas,.1:.1:.9,.1);
    PlotScalpMaps(sm,chanlocs,[],.1:.1:.9,suffixes_in);
end