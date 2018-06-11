% ClassifyWithPupilSize_v3
%
% Created 7/11/13 by DJ for one-time use.
% Updated 4/29/16 by DJ - added code from GetTagScoresForEachObject.m to
%  get each object's TAG score and save it

%% Set up
homedir = '/Users/jangrawdc/Documents/LiincData/3DSearch';
folders = {'2013-03-27-Pilot-S22', '2013-03-28-Pilot-S23', ...
    '2013-04-29-Pilot-S24', '2013-05-01-Pilot-S25', '2013-05-02-Pilot-S26',...
    '2013-05-03-Pilot-S27', '2013-05-07-Pilot-S28', '2013-05-10-Pilot-S29',...
    '2013-05-15-Pilot-S30', '2013-06-06-Pilot-S32'};
subjects = 22;%[22:30 32];
sessions_cell = {2:14, [3 6:17], 1:15, 1:15, 1:15, 1:15, 1:15, 1:15, [1:10 12:15], 2:16};
offsets = [-12 -4 48 60 68 68 92 112 -32 88];
epochRange = [-1000 4000];
cvmode = '10fold';
cd(homedir);
load('TEMP_dwellFeature');

%%%%% PICK FEATURES %%%%%
useEEG = 1;
usePS = 1;
useDT = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PICK NTRIALS %%%%%
nSessionsToInclude = nan; % NaN for all sessions
%%%%%%%%%%%%%%%%%%%%%%%%%

% Declare PS bins
tBaseline = [-1000 0];%[0 100];
binwidth = 500; % in ms
binstart = 0:500:2500; % in ms

% Declare loop variables
xeye_all = cell(1,numel(subjects));
ps_all = cell(1,numel(subjects));
Az_ps_all = [];
ps_median_targ_all = [];
ps_median_dist_all = [];
clear R Az

for iSubj = 1:numel(subjects);
    subject = subjects(iSubj);
    sessions = sessions_cell{iSubj};
    % crop sessions if specified
    if ~isnan(nSessionsToInclude) 
        sessions = sessions(1:nSessionsToInclude);
    end
    offset = offsets(iSubj);
    
    %% Load EEG and behavior structs
    y = loadBehaviorData(subject,sessions,'3DS');
    EEG = pop_loadset('filename',...
        sprintf('3DS-%d-all-filtered-noduds-noeog-epoched-ica.set',subject));
    % crop epochs if specified
    if ~isnan(nSessionsToInclude) 
        isNewEpoch = [true, diff([EEG.event.epoch])>0];
        isNewSession = diff([EEG.event(isNewEpoch).init_index])<0;        
        epochSession = [1, cumsum(isNewSession)+1];
        trialrej = epochSession>nSessionsToInclude;
        EEG = pop_rejepoch(EEG,trialrej);
    end
    %% Get eye pos and pupil size
    [ps,xeye,yeye] = deal(cell(1,numel(sessions)));
    for i=1:numel(sessions)        
        load(sprintf('3DS-%d-%d-eyepos',subject,sessions(i)));
        xeye{i} = eyepos(:,1);
        yeye{i} = eyepos(:,2);
        ps{i} = InterpolateBlinks(pupilsize,y(i).eyelink.record_time-1+(1:length(pupilsize)),y(i));
    end 
    
    %% Get epochs
    [ps_epoch,t_epoch,isTargetEpoch] = GetEpochedPupilSize(y,ps,EEG,epochRange);
    xeye_epoch = GetEpochedPupilSize(y,xeye,EEG,epochRange);
    yeye_epoch = GetEpochedPupilSize(y,yeye,EEG,epochRange);

    %% Normalize pupil size to be percentage change
    mean_ps = nanmean(cat(1,ps{:}));
    ps_epoch_pct = ps_epoch/mean_ps*100;

    %% Remove baseline
    isBase = t_epoch>=tBaseline(1) & t_epoch<=tBaseline(2);
    baseline = nanmean(ps_epoch_pct(:,isBase),2);
    ps_epoch_pct_baserem = ps_epoch_pct - repmat(baseline,1,size(ps_epoch,2));

    %% Get AUC values
    % Build windowed classifier
    binstart_samples = find(ismember(t_epoch,binstart));    
    binwidth_samples = round(binwidth/(t_epoch(2)-t_epoch(1)));
    bindata = nan(size(ps_epoch_pct_baserem,1),numel(binstart));
    for i=1:numel(binstart)        
        bindata(:,i) = nanmean(ps_epoch_pct_baserem(:,binstart_samples(i)+(1:binwidth_samples)-1),2);        
    end
    bindata(isnan(bindata)) = 0; % if unknown, say no info.
    load(sprintf('ALLEEG%d_NoeogEpochedIcaCropped.mat',subject)); % ALLEEG
    % crop epochs if specified
    if ~isnan(nSessionsToInclude) 
        for k=1:2
            isNewEpoch = [true, diff([ALLEEG(k).event.epoch])>0];
            isNewSession = diff([ALLEEG(k).event(isNewEpoch).init_index])<0;        
            epochSession = [1, cumsum(isNewSession)+1];
            trialrej = epochSession>nSessionsToInclude;
            ALLEEG(k) = pop_rejepoch(ALLEEG(k),trialrej);
        end
    end
    
    % Add in de-meaned dwell time
    dwellColumn = dwellFeature{iSubj}; % de-meaned dwell time
    combodata = [];
    if usePS
        combodata = [combodata, bindata];
    end
    if useDT
        combodata = [combodata, dwellColumn];
    end
    % Classify!
    R(iSubj) = ClassifyWithEegAndDwellTime(ALLEEG,y,cvmode,offset,combodata,useEEG);
    
    Az(iSubj) = R(iSubj).Az;

    %% Get ps plot info

    % Get medians & stderrors    
    ps_median_targ = nanmedian(ps_epoch_pct_baserem(isTargetEpoch==1,:),1);
    ps_stderr_targ = nanstd(ps_epoch_pct_baserem(isTargetEpoch==1,:),[],1)/sqrt(sum(isTargetEpoch==1));
    ps_median_dist = nanmedian(ps_epoch_pct_baserem(isTargetEpoch==0,:),1);
    ps_stderr_dist = nanstd(ps_epoch_pct_baserem(isTargetEpoch==0,:),[],1)/sqrt(sum(isTargetEpoch==0));

    % Add to results
    ps_median_targ_all(iSubj,:) = ps_median_targ;
    ps_median_dist_all(iSubj,:) = ps_median_dist;
    ps_all{iSubj} = ps_epoch_pct_baserem;
    xeye_all{iSubj} = xeye_epoch;
    yeye_all{iSubj} = yeye_epoch;
end

%% Get predicted target objects
iObjects_hci_pt = cell(1,numel(subjects));
for i=1:numel(subjects)
    fprintf('--- Subject %d/%d ---\n',i,numel(subjects));
    load(sprintf('ALLEEG%d_NoeogEpochedIcaCropped.mat',subjects(i))); % ALLEEG
    iObjects_hci_pt{i} = GetEegPredictedTargets(subjects(i),sessions_cell{i},ALLEEG,R(i).y,R(i).truth);
end

%% Get TSP and stats
clear stats
usegridconstraints = true;
levelname = 'GridHuge.jpg';
fracSensitivity = 0.9999999; % used to calculate nSensitivity
for i=1:numel(subjects)
    fprintf('--- Subject %d/%d ---\n',i,numel(subjects));
    subject = subjects(i);    
    sessions = sessions_cell{i};
    [stats(i),TagTours{i},plotInfos{i}] = CalculateImprovement(subject,sessions,levelname,iObjects_hci_pt{i},usegridconstraints,fracSensitivity);    
end

%% SET CHANCE LEVELS
chance.precision = [25 25];
chance.pctFound = [20/56*.25*100 25];
chance.pctDistance = 100;




%% %% Get object lists
[objName, objTagScore, objLocation, targetCat] = deal(cell(1,numel(subjects)));
eyelinkoreeg='eeg';
nRows = ceil(sqrt(numel(subjects)));
nCols = ceil(numel(subjects)/nRows);
xHist = linspace(0,0.07,100);

for iSubject=1:numel(subjects)
    cd(fullfile(homedir,folders{iSubject}));
    % find all object names and categories
    subject = subjects(iSubject);
    sessions = sessions_cell{iSubject};
    [objects, objectnames, objectlocations, objecttimes, objectisessions] = GetObjectList(subject,sessions,eyelinkoreeg);    
    % find object categories    
    iCat = nan(1,numel(objectnames));
    for i=1:numel(categories)
        isInCat = strncmp(categories{i},objectnames,length(categories{i}));
        iCat(isInCat) = i;
    end
    
    objName{iSubject} = objectnames';
    objTagScore{iSubject} = plotInfos{iSubject}.tagScore';
    objLocation{iSubject} = plotInfos{iSubject}.objlocs;
    
    % Infer target category
    iTargetCat = unique(iCat(plotInfos{iSubject}.iTargets));
    targetCat{iSubject} = categories{iTargetCat};

    % Make histogram for each subject/category
    xHist = linspace(0,max(objTagScore{iSubject}),50);
    histCatTagScore = zeros(numel(xHist),numel(categories));
    legendstr = cell(size(categories));
    for i=1:numel(categories)
        histCatTagScore(:,i) = hist(objTagScore{iSubject}(iCat==i),xHist);
        if i==iTargetCat
            legendstr{i} = [categories{i} '*'];
        else
            legendstr{i} = categories{i};
        end
    end

    % Plot results
    subplot(nRows,nCols,iSubject);
    plot(histCatTagScore);
    legend(legendstr,'interpreter','none');
    xlabel('TAG score')
    ylabel('# objects')
    title(sprintf('Subject %d histogram',subject));
    
end




%% Save results
README = 'objName{i}{j} is the name of object #j in the environment of subject number subjects(i). objLocation{i}(j,:) is the x,y location of the object for plotting. objTagScore{i}(j) is the TAG score for object number j from the classifier of subject number subjects(i). The instructed target category for that subject was targetCat{i}.';
suffix = '';
if useEEG
    suffix = [suffix 'eeg'];
end
if usePS
    suffix = [suffix 'ps'];
end
if useDT
    suffix = [suffix 'dt'];
end
if isnan(nSessionsToInclude)
    fprintf('Saving %s/TagScoresForEachObject_%s.mat...\n',homedir,suffix);
    save(sprintf('%s/TagScoresForEachObject_%s.mat',homedir,suffix),'subjects','categories','objName','objTagScore','targetCat','README');
else
    fprintf('Saving %s/TagScoresForEachObject_%s_%02dsess.mat...\n',homedir,suffix,nSessionsToInclude);
    save(sprintf('%s/TagScoresForEachObject_%s_%02dsess.mat',homedir,suffix,nSessionsToInclude),'subjects','categories','objName','objTagScore','targetCat','README');
end
fprintf('Done!\n');
