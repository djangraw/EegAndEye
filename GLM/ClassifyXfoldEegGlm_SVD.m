function [prediction, truth] = ClassifyXfoldEegGlm_SVD(subject,nFolds,nEventsPerTrial,tRange,components)

% Attempts to classify each trial in a squares experiment using SVD results
%
% [prediction, truth] = ClassifyXfoldEegGlm_SVD(subject,nFolds)
%
% INPUTS:
% - subject is a scalar indicating the subject number.
% - nFolds is a scalar indicating the number of structs you would like to
% create.
% - nEventsPerTrial is a scalar indicating the number of saccade events per
% trial (default = 5).
% - tRange is a 2-element vector indicating the start and end times (in ms) 
% that should be used for the SVD.
% - components is a vector of the SVD components that should be used for
% classification.
%
% OUTPUTS:
% - prediction is the predicted sequence for each trial (1 row = 1 trial)
% - truth is the actual sequence for each trisl
%
% Is careful to recommend separations between trials.  Usually called by
% RunXfoldEegGlm().
%
% Created 5/15/12 by DJ based on ClassifyXfoldEegGlm.m.
% Updated 3/22/13 by DJ - use asymmetric influence, GetGlmRegressors_v2p0
% Updated 4/30/13 by DJ - added support for responseFns in cells

if nargin<3 || isempty(nEventsPerTrial)
    nEventsPerTrial = 5;
end
if nargin<4 || isempty(tRange)
    tRange = [200 500];
end
if nargin<5 || isempty(components)
    components = 1;
end

% Options
do_plots = 0;
do_lines = 0;
chanstoplot = 1:numel(components);

% Initialize variables
prediction = cell(1,nFolds);
truth = cell(1,nFolds);
accuracy = zeros(1,nFolds);
seq_accuracy = zeros(1,nFolds);
nTrials = zeros(1,nFolds);

saccade_events = {'Dist-0T','Dist-1T','Integ','Compl','Irrel','Circle-D','Circle-T'};

for iFold=1:nFolds
    % Load
    TRAIN = load(sprintf('sq-%d-GLMresults-Nuisance-Training-fold%d.mat',subject,iFold));
    TEST = load(sprintf('sq-%d-EEGtest-Nuisance-Testing-fold%d.mat',subject,iFold));
    % Unpack
%     nuisance_events = TRAIN.nuisance_events;
    if ~iscell(TRAIN.regressor_events{1})       
        regressor_events = TRAIN.regressor_events;
    else
        regressor_events = TRAIN.regressor_events{TRAIN.iLevel};
    end
    event_influence_ms = TRAIN.influence; % [min max] in ms
    if isnumeric(TRAIN.responseFns)
        responseFns = TRAIN.responseFns;
    else
        responseFns = TRAIN.responseFns{TRAIN.iLevel};
    end
    tResponse = TRAIN.tResponse;
    EEG = TEST.EEGtest;    
    
    % Run SVD
    weights = ApplySvdToGlmResults(responseFns,tResponse,tRange,EEG.chanlocs,regressor_events);
    % Get timecourses
    timecourse = zeros(length(components),length(tResponse),length(regressor_events));
    for i=1:numel(regressor_events)
        for j=1:numel(components)
            timecourse(j,:,i) = weights(:,components(j))'*responseFns(:,:,i); % use dominant component weights(:,1)
        end
    end
    
    % Find constants
    dt = 1000/EEG.srate;
    regressor_range = round(event_influence_ms/dt); % how many samples should each response function extend?
    % Handle 1-element ranges
    if numel(regressor_range) == 1
        regressor_range = [-regressor_range, regressor_range];
    end
    t = (1:EEG.pnts)*dt; % time vector for EEG

    % Find event times
    Nr = numel(regressor_events);
    event_times = cell(1,Nr);
    for i=1:Nr
        if ismember(regressor_events{i},saccade_events)
            event_times{i} = [EEG.event(strcmp(regressor_events{i},{EEG.event(:).type})).latency]*dt;
        else % exclude button events
            event_times{i} = [];
        end
    end

    % Find blink events
%     artifact_times = [EEG.event(ismember({EEG.event(:).type}, artifact_events)).latency]*dt;

    % Get trial times
    tTrialStart = [EEG.event(strncmp('TrialStart',{EEG.event.type},length('TrialStart'))).latency]*dt;
        
    disp('Getting regressors...');
    % Get regressors
%     s = GetGlmRegressors(t,event_times,[],regressor_range);
    s = GetGlmRegressors_v2p0(t,event_times,[],regressor_range);
%     sClean = GetGlmRegressors(t,event_times,artifact_times,regressor_range);

    % Get event times for each trial
    iEvents = zeros(numel(tTrialStart),nEventsPerTrial); 
    eventTypes = zeros(numel(tTrialStart),nEventsPerTrial); 
    for iTrial = 1:numel(tTrialStart)-1
        iThis = find(sum(s,1) & t>tTrialStart(iTrial) & t<tTrialStart(iTrial+1));
        if numel(iThis)==nEventsPerTrial
            iEvents(iTrial,:) = iThis;
            [eventTypes(iTrial,:),~] = find(s(:,iThis));
        end
    end
    iThis = find(sum(s,1) & t>tTrialStart(end));
    if numel(iThis)==nEventsPerTrial
        iEvents(end,:) = iThis;
    end
    % Crop to only include complete trials
    isComplete = sum(iEvents,2)>0;
    iEvents = iEvents(isComplete,:);
    eventTypes = eventTypes(isComplete,:);
    tTrialStart = tTrialStart(isComplete);
    
    % Get possibe sequences
    [sequences, chance] = GetAllTrialSequences(regressor_events); %#ok<NASGU> (suppress compiler warning)
    % Compare actual data to sequence templates
    iSeq = zeros(1,numel(tTrialStart));
    iTrue = zeros(1,numel(tTrialStart));
    seqError = zeros(numel(tTrialStart),size(sequences,1));
    fitError = zeros(1,numel(tTrialStart));
    for i = 1:numel(tTrialStart)-1 % trial
        iTrialStart = iEvents(i,1)+regressor_range(1)-1;
        iTrialEnd = iEvents(i,end)+regressor_range(2)+1;
        data = weights(:,components)'*EEG.data(:,iTrialStart:iTrialEnd);        
        % Plot true eeg
        if do_plots
            figure(1); clf; hold on; 
            imagesc(data(chansToPlot,:)); PlotVerticalLines(iEvents(i,:)-iTrialStart,'k--');        
            set(gca,'clim',[-30 30]); colorbar; 
            title(['Fold/Trial ' num2str(iFold) '/' num2str(i) ' TRUE sequence = ' num2str(eventTypes(i,:))]);
            figure(2);
        end
        if do_lines
            figure(3); clf; 
            for k=1:numel(chanstoplot)
                subplot(numel(chanstoplot),1,k); hold on; ylabel(sprintf('Component %d activity', components(chanstoplot(k))));
                plot((1:size(data,2))*dt, data(chanstoplot(k),:)/sqrt(var(data(chanstoplot(k),:))),'b'); 
            end
%             subplot(3,1,1); hold on; ylabel(EEG.chanlocs(17).labels);          
%             plot((1:size(data,2))*dt, data(17,:)/sqrt(var(data(17,:))),'b'); 
%             subplot(3,1,2); hold on; ylabel(EEG.chanlocs(38).labels);
%             plot((1:size(data,2))*dt, data(38,:)/sqrt(var(data(38,:))),'b'); 
%             subplot(3,1,3); hold on; ylabel(EEG.chanlocs(58).labels);
%             plot((1:size(data,2))*dt, data(58,:)/sqrt(var(data(58,:))),'b'); 
        end
        iTrue(i) = find(ismember(sequences,eventTypes(i,:), 'rows'),1);        
        for j = 1:size(sequences,1) % sequence
            template = zeros(size(data));
            for k = 1:nEventsPerTrial % saccade
                template(:,(iEvents(i,k)-iTrialStart+regressor_range(1)):(iEvents(i,k)-iTrialStart+regressor_range(2))) = ...
                    template(:,(iEvents(i,k)-iTrialStart+regressor_range(1)):(iEvents(i,k)-iTrialStart++regressor_range(2))) + ...
                    timecourse(:,:,sequences(j,k));
            end
            template(template==0) = NaN;            
            if do_lines && j==iTrue(i) 
                for k=1:numel(chanstoplot)
                    subplot(numel(chanstoplot),1,k); hold on;
                    plot((1:size(data,2))*dt, template(chanstoplot(k),:)/sqrt(nanvar(template(chanstoplot(k),:))),'c');
                end
%                 subplot(3,1,1); hold on;                
%                 plot((1:size(data,2))*dt, template(17,:)/sqrt(nanvar(template(17,:))),'r');
%                 subplot(3,1,2); hold on; 
%                 plot((1:size(data,2))*dt, template(38,:)/sqrt(nanvar(template(38,:))),'r');
%                 subplot(3,1,3); hold on; 
%                 plot((1:size(data,2))*dt, template(58,:)/sqrt(nanvar(template(58,:))),'r');
            end
            seqError(i,j) = nanmean(template(:).*data(:));
            if do_plots
                subplot(4,4,j); cla; hold on;
                imagesc(template(chansToPlot,:)); PlotVerticalLines(iEvents(i,:)-iTrialStart,'k--');
                set(gca,'clim',[-10 10]); colorbar; title(['sequence = ' num2str(sequences(j,:)) ', error = ' num2str(seqError(i,j))]);            
            end
        end        
        [fitError(i),iSeq(i)] = max(seqError(i,:));
        if do_plots
            % plot labels
            subplot(4,4,iSeq(i)); text(0,3,' *WINNER*')
            subplot(4,4,iTrue(i)); text(0,1,' *TRUE*')
        end
        if do_lines
            template = zeros(size(data));
            for k = 1:nEventsPerTrial % saccade
                template(:,(iEvents(i,k)-iTrialStart+regressor_range(1)):(iEvents(i,k)-iTrialStart+regressor_range(2))) = ...
                    template(:,(iEvents(i,k)-iTrialStart+regressor_range(1)):(iEvents(i,k)-iTrialStart+regressor_range(2))) + ...
                    timecourse(:,:,sequences(iSeq(i),k));
            end
            template(template==0) = NaN;
            subplot(numel(chanstoplot),1,1); hold on;
            title(sprintf('Fold/Trial %d/%d,  TRUE sequence = %s\n PREDICTED sequence = %s',...
                iFold,i,num2str(eventTypes(i,:)),num2str(sequences(iSeq(i),:))));
%             plot((1:size(data,2))*dt, template(17,:)/sqrt(nanvar(template(17,:))),'g'); 
            for k=1:numel(chanstoplot)
                subplot(numel(chanstoplot),1,k)
                plot((1:size(data,2))*dt, template(chanstoplot(k),:)/sqrt(nanvar(template(chanstoplot(k),:))),'r');
                PlotVerticalLines((iEvents(i,:)-iTrialStart)*dt,'k--');            
            end
%             PlotVerticalLines((iEvents(i,:)-iTrialStart)*dt,'k--');            
%             subplot(3,1,2); hold on;
% %             plot((1:size(data,2))*dt, template(38,:)/sqrt(nanvar(template(38,:))),'g'); 
%             PlotVerticalLines((iEvents(i,:)-iTrialStart)*dt,'k--');
%             subplot(3,1,3); hold on;
% %             plot((1:size(data,2))*dt, template(58,:)/sqrt(nanvar(template(58,:))),'g'); 
%             PlotVerticalLines((iEvents(i,:)-iTrialStart)*dt,'k--');         
            legend('data','template','predicted','events');
            
        end
        if do_plots || do_lines
            % get input on whether to keep going
            foo = input('','s');
            if ~isempty(foo) % stop if user enters text
                return;
            end
        end
        
    end
    % Create metrics of success
    prediction{iFold} = sequences(iSeq(1:end-1),:);
    truth{iFold} = eventTypes(1:end-1,:);
    iCompl = find(strcmp('Compl',regressor_events));
    accuracy(iFold) = mean(any(prediction{iFold}==iCompl,2)==any(truth{iFold}==iCompl,2));
    seq_accuracy(iFold) = mean(iSeq==iTrue);
    nTrials(iFold) = size(truth{iFold},1);
    fprintf('class accuracy: %d/%d = %f\n',round(accuracy(iFold)*nTrials(iFold)), nTrials(iFold), accuracy(iFold));
    fprintf('sequence accuracy: %d/%d = %f\n',round(seq_accuracy(iFold)*nTrials(iFold)), nTrials(iFold), seq_accuracy(iFold));
end
% Display total accuracy
totalAccuracy = accuracy*nTrials'/sum(nTrials);
totalSeqAccuracy = seq_accuracy*nTrials'/sum(nTrials);
fprintf('TOTAL class accuracy: %d/%d = %f\n',round(accuracy*nTrials'), sum(nTrials), totalAccuracy);
fprintf('TOTAL sequence accuracy: %d/%d = %f\n',round(seq_accuracy*nTrials'), sum(nTrials), totalSeqAccuracy);