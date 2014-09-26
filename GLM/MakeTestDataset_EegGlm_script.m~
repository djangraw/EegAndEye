% MakeTestDataset_EegGlm_script.m
%
% Created 4/22/13 by DJ.
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 3/19/14 by DJ - tResponse, influence, artifact_influence in cells

%% Set up
N = 3; % subjects
D = 1; % electrodes
T = 20; % offsets
M = 2; % events
eventnames = cell(1,M);
for i=1:M
    eventnames{i} = sprintf('event%d',i);
end
clear chanlocs
for i=1:D
    chanlocs(i).labels = sprintf('chan%d',i);
end

Fs = 100; % sampling rate in Hz
n = (1:N)*10000; % time points for each subject
trueBg = [(1:T)-T/2; repmat(T/2,1,T/2), T/2:-1:1];
% trueBg = reshape(1:(M*T),T,M)'; % size of each event's effect (group level)
% trueBg = trueBg - mean(trueBg(:)); % demean
trueCg = trueBg(2,:)-trueBg(1,:);
figure(111); plot((0:T-1)/Fs*1000,trueCg);
xlabel('time (ms)')
ylabel('Contrast')
title('True group contrast')
% make true beta matrix
trueB = zeros(M,T,N);
betaSize = [3 2 1];
betaSize = betaSize/mean(betaSize); % normalize
for i=1:N
    trueB(:,:,i) = trueBg*betaSize(i);    
end
trueSigma = 20*([1 1 1]); % variance of each effect

%% create
eventtimes = cell(1,N);
eventtypes = cell(1,N);
clear ALLEEG;
for i=1:N    
    % make events
    eventtimes = (1:T:(n(i)-T))'; % (in samples) regularly, with no breaks
    eventtypes = ceil(rand(size(eventtimes))*2);
    % make data
    eegdata = zeros(D,n(i));    
    for j=1:numel(eventtimes)
        eegdata(:,eventtimes(j)+(0:T-1)) = eegdata(:,eventtimes(j)+(0:T-1)) + trueB(eventtypes(j),:,i);
    end        
    % add random noise
    eegdata = eegdata+randn(size(eegdata))*trueSigma(i);    
    
    % Set up event matrix
    eventmat = [eventnames(eventtypes)',mat2cell(eventtimes/Fs, ones(length(eventtimes),1), 1)];
    eventmat = [eventmat; {'epochevent', 1/Fs; 'epochevent', n(i)/2/Fs}];
    % Add to EEGLAB struct
    EEG = pop_importdata('data','eegdata','srate',Fs,'setname',sprintf('test-%d',i),'subject',i);
    EEG = pop_importevent( EEG, 'event','eventmat','fields',{'type','latency'},'append','no','timeunit',1);
    EEG = eeg_checkset(EEG); 
    % Make epochs
    EEG = pop_epoch(EEG, {'epochevent'},[EEG.xmin EEG.xmax/2]);
    % Add weird fields
    EEG.etc.ureventweights = ones(size(EEG.urevent));
    EEG.chanlocs = chanlocs;
    % Save results
    pop_saveset(EEG,'filename',EEG.setname);
    ALLEEG(i) = EEG;
end
disp('---Done!')

%% Solve 1st-level
% Set up
regressor_events = {eventnames};
influence = {[0 (T-1)/Fs*1000]}; % in ms
artifact_events = {};
artifact_influence = {[]};
offset = 0;
stddev = 0;
vthresh = inf;
method = 'mvregress';
trial_rej_rules = {};
iLevel = 1;

% Run
for i=1:N        
    fprintf('Subject %d...\n',i);
    EEG = ALLEEG(i);
    filenames = {sprintf('test-%d-GLMresults-test',i)};
    dataset = EEG.setname;
    
    % Solve
    [newRF, newTR] = RunEegGlm(EEG,regressor_events{iLevel},influence{iLevel},artifact_events,artifact_influence{iLevel},method,stddev);
    responseFns{iLevel} = newRF;
    tResponse{iLevel} = newTR;
    
    if ~isempty(filenames{iLevel})
        % Update status
        fprintf('%s Saving Results as %s...\n',datestr(now,16), filenames{iLevel});
        % Save results
        save(filenames{iLevel},'EEG','responseFns','tResponse','regressor_events',...
            'filenames','iLevel','artifact_events','dataset','offset','influence',...
            'artifact_influence','stddev','vthresh','method','trial_rej_rules');        
    end
    disp('Done!')
%     EEG = SubtractOutGlmResponses(EEG,responseFns,influence,regressor_events{iLevel},stddev);

end

%% Perform 2nd-level GLM
% Set up
prefix = 'test';
subjects = 1:N;
glmType = 'test';
% Set up contrasts
plus_events = eventnames{2};
minus_events = eventnames{1};
baseline_win = [0 50];
multcorrect = 'fdr';

% Load results
for i=1:numel(subjects)
    fprintf('Loading subject %d/%d...\n',i,numel(subjects));
    filename = sprintf('%s-%d-GLMresults-%s',prefix,subjects(i),glmType);
    R(i) = load(filename);
end

% Set up betas
[contrastFns, contrastVar] = SetUpTopLevelGlm(R,plus_events,minus_events,baseline_win);

% Run top-level GLM
[group_RF,group_P] = RunTopLevelGlm_EEG(contrastFns,contrastVar,multcorrect);
group_Z = norminv(group_P);

for i=1:N
    figure(100+i); subplot(3,1,1); set(gca,'ylimmode','auto');
end
disp('---Done!')

%% Plot results
if iscell(R(1).tResponse)
    tResponse = R(1).tResponse{R(1).iLevel};
else
    tResponse = R(1).tResponse;
end
figure(208); clf;
MakeFigureTitle(sprintf('%s GLM, %s - %s contrast',glmType,plus_events,minus_events));

subplot(2,2,1); hold on;
plot(tResponse, mean(contrastFns,4));
plot(tResponse, mean(mean(contrastFns,4),1),'k','linewidth',2);
xlabel('time (ms)')
ylabel('voltage (uV)')
title('Subject Average')

subplot(2,2,2); hold on;
plot(tResponse, group_RF);
plot(tResponse, mean(group_RF,1),'k','linewidth',2);
xlabel('time (ms)')
ylabel('voltage (uV)')
title('Hierarchical GLM result')

subplot(2,2,3); hold on;
plot(tResponse, mean(contrastFns,4)-group_RF);
plot(tResponse, mean(mean(contrastFns,4)-group_RF,1),'k','linewidth',2);
xlabel('time (ms)')
ylabel('voltage (uV)')
title('Difference')

subplot(2,2,4); hold on;
plot(tResponse, group_P);
plot(tResponse, mean(group_P,1),'k','linewidth',2);
isSignifHi = any(group_P>0.975,1);
isSignifLo = any(group_P<0.025,1);
plot(get(gca,'xlim'),[0 0],'k--')
plot(get(gca,'xlim'),[1 1],'k--')
plot(get(gca,'xlim'),[0.025 0.025],'k:')
plot(get(gca,'xlim'),[.975 .975],'k:')
plot(tResponse(isSignifHi),repmat(1.1,1,sum(isSignifHi)),'g.');
plot(tResponse(isSignifLo),repmat(-0.1,1,sum(isSignifLo)),'g.');
ylim([-0.1 1.1]);
xlabel('time (ms)')
ylabel('P value')
title(sprintf('Hierarchical GLM p values (mc correction: %s)',multcorrect))
