% SynthDataGlmTest_script
%
% Created 3/24/14 by DJ.

% Declare parameters
nTrials = 10;
sqPerTrial = 3;
eventList = {'event1';'event2';'event3'};
Fs = 100;
D = 4;
h = repmat(cat(3,1:10, 10:-1:1, [ones(1,5),-ones(1,5)]),[D 1 1]);
T = size(h,2);

nSubjects = 3;
A_subject = 10*(1:nSubjects);
noise_std = 1;

% Declare chanlocs
load SampleChanlocs
chanlabels = {'FZ','CZ','PZ','OZ'};
chanlocs = chanlocs(ismember({chanlocs.labels},chanlabels));

for iSubj=1:nSubjects
    % Create events list
    eventTypes = repmat(eventList,10,1);
    eventTimes = zeros(sqPerTrial*nTrials,1);
    for i=1:nTrials
        eventTimes((i-1)*sqPerTrial + (1:sqPerTrial)) = 5*i + (1:sqPerTrial);
    end

    % Construct data
    data = zeros(D,Fs*(eventTimes(end)+5));
    t = (1:size(data,2))/Fs; % time in s
    for i=1:numel(eventTimes)
        iTime = find(t>eventTimes(i),1);
        iEventType = find(strcmp(eventTypes{i},eventList));        
        data(:,iTime-1+(1:T)) = data(:,iTime-1+(1:T)) + A_subject(iSubj)*h(:,:,iEventType);
%         data(:,iTime-1+(1:T)) = data(:,iTime-1+(1:T)) + h(:,:,iEventType)* A_subject(iSubj)*noise_std .* randn(size(h(:,:,iEventType))); % constant SNR
    end
    % Add noise
    data = data + A_subject(iSubj)*noise_std*randn(size(data)); % constant noise level

    % Import into EEGLAB struct
    events = [eventTypes, num2cell(eventTimes)];
    EEG = pop_importdata('data',data,'srate',Fs);
    EEG.chanlocs = chanlocs;
    EEG = eeg_checkset(EEG);
    EEG = pop_importevent(EEG,'event','events','fields',{'type','latency'});
    EEG = pop_epoch(EEG,{'event1'},[-1 4]);


    % Run Level 1 GLM
    EEG.etc.ureventweights = ones(numel(EEG.event));

    regressor_events = {eventList};
    influence = {[0 T-1]/Fs*1000};% in ms
    artifact_events = {''};
    artifact_influence = {[0 -1]};
    method = 'mvregress';
    offset = 0;
    stddev = 0;
    dataset = sprintf('TEST-%d',iSubj);
    vthresh = inf;
    trial_rej_rules = {''};
    filenames = {sprintf('TEST-%d.mat',iSubj)};

    nLevels = length(regressor_events);
    [responseFns, tResponse] = deal(cell(1,nLevels));
    EEGnew = EEG; % create copy to be modified during run (but only original will be saved)
    % Run analysis
    for iLevel=1:nLevels % for each level of analysis
        % Update status
        fprintf('%s - Running Analysis Level %d/%d...\n',datestr(now,16),iLevel,nLevels);

        [newRF, newTR] = RunEegGlm(EEGnew,regressor_events{iLevel},influence{iLevel},artifact_events,artifact_influence{iLevel},method,stddev);
        responseFns{iLevel} = newRF;
        tResponse{iLevel} = newTR;
    %     responseFns = []; tResponse  = [];

        if ~isempty(filenames{iLevel})
            % Update status
            fprintf('%s Saving Results as %s...\n',datestr(now,16), filenames{iLevel});        
            % Save results
            save(filenames{iLevel},'EEG','responseFns','tResponse','regressor_events',...
                'filenames','iLevel','artifact_events','dataset','offset','influence',...
                'artifact_influence','stddev','vthresh','method','trial_rej_rules');
        end
        EEGnew = SubtractOutGlmResponses(EEGnew,responseFns{iLevel},influence{iLevel},regressor_events{iLevel},stddev);
    end

end
    
% Plot results
% 
% figure;
% nEvents = size(h,3);
% for i=1:nEvents
%     subplot(nEvents,2,2*i-1);
%     imagesc(h(:,:,i));
%     colorbar
%     subplot(nEvents,2,2*i);
%     imagesc(responseFns{1}(:,:,i));
%     colorbar
% end

%% Top-level GLM

eventWeights = eye(nEvents);

% Load results
for iSubj = 1:nSubjects
    R(iSubj) = load(sprintf('TEST-%d.mat',iSubj));
end

% Zscores
[contrastFns, contrastVar, contrastZ]  = SetUpTopLevelGlm_flex(R,eventList,eventWeights);
%
multcorrect = 'none';
% run level 2
[group_RF,group_P] = RunTopLevelGlm_EEG(contrastFns,contrastVar,multcorrect);

%% Plot results
group_Z = norminv(group_P);

figure;
nEvents = size(h,3);
for i=1:nEvents
    subplot(nEvents,3,3*i-2);
    imagesc(h(:,:,i)*mean(A_subject));
    colorbar
    subplot(nEvents,3,3*i-1);
    imagesc(group_RF(:,:,i));
    colorbar
    subplot(nEvents,3,3*i);
    imagesc(group_Z(:,:,i));
    colorbar
end

