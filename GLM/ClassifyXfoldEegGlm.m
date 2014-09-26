function [prediction, truth] = ClassifyXfoldEegGlm(subject,nFolds,nEventsPerTrial)

% Attempts to classify each trial in a squares experiment.
%
% [prediction, truth] = ClassifyXfoldEegGlm(subject,nFolds,nEventsPerTrial)
%
% INPUTS:
% - subject is a scalar indicating the subject number.
% - nFolds is a scalar indicating the number of structs you would like to
% create.
% - nEventsPerTrial is a scalar indicating the number of saccade events per
% trial (default = 5).
%
% OUTPUTS:
% - prediction is the predicted sequence for each trial (1 row = 1 trial)
% - truth is the actual sequence for each trisl
%
% Is careful to recommend separations between trials.  Usually called by
% RunXfoldEegGlm().
%
% Created 2/1/12 by DJ.
% Updated 3/19/12 by DJ - comments
% Updated 5/15/12 by DJ - added nEventsPerTrial input, v1pt5 support
% Updated 5/16/12 by DJ - added match_method, use_posterior,
% pooled_complete parameters, bar plot as part of do_lines
% Updated 3/22/13 by DJ - use asymmetric influence, GetGlmRegressors_v2p0
% Updated 4/30/13 by DJ - added support for responseFns in cells

% Handle inputs
if nargin<3
    nEventsPerTrial = 5;
end

% Options
do_plots = 0;
do_lines = 0;
chanstoplot = [17, 38, 58];
match_method = 'AvgRefNormProduct'; % Correlation methods: 'Product' 'NormProduct' 'Difference' 'Pearson' 'Spearman'
use_posterior = 1; % multiply match score by chance of each sequence occurring
pooled_complete = 0; % use the pooled probability of all complete trials to predict completeness

% Initialize variables
prediction = cell(1,nFolds);
truth = cell(1,nFolds);
accuracy = zeros(1,nFolds);
seq_accuracy = zeros(1,nFolds);
nTrials = zeros(1,nFolds);

saccade_events = {'Dist-0T','Dist-1T','Integ','Compl','Irrel'};

for iFold=1:nFolds
    fprintf('Fold %d/%d...\n',iFold,nFolds);
    
    % Load
    TRAIN = load(sprintf('sq-%d-GLMresults-Nuisance-Training-fold%d.mat',subject,iFold));
    TEST = load(sprintf('sq-%d-EEGtest-Nuisance-Testing-fold%d.mat',subject,iFold));
    
    % Unpack
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

    disp('Getting regressors...');
    % Get regressors
%     s = GetGlmRegressors(t,event_times,[],regressor_range);
    s = GetGlmRegressors_v2p0(t,event_times,[],regressor_range);

%     % Get regressors with artifacts removed
%     artifact_times = [EEG.event(ismember({EEG.event(:).type}, artifact_events)).latency]*dt;
%     sClean = GetGlmRegressors(t,event_times,artifact_times,Nt);

    % Get event times for each trial
    tTrialStart = [EEG.event(strncmp('TrialStart',{EEG.event.type},length('TrialStart'))).latency]*dt; % Get trial times    
    nTrials(iFold) = numel(tTrialStart);
    iEvents = zeros(nTrials(iFold),nEventsPerTrial); 
    eventTypes = zeros(nTrials(iFold),nEventsPerTrial); 
    for iTrial = 1:nTrials(iFold) % for each trial
        if iTrial==nTrials(iFold)
            iThis = find(sum(s,1) & t>tTrialStart(end));
        else
            iThis = find(sum(s,1) & t>tTrialStart(iTrial) & t<tTrialStart(iTrial+1));
        end
        if numel(iThis)==nEventsPerTrial
            iEvents(iTrial,:) = iThis;
            [eventTypes(iTrial,:),~] = find(s(:,iThis));
        end
    end    
    
    % Crop to only include well-done (viewed all 5 squares) trials
    isWellDone = sum(iEvents,2)>0;
    iEvents = iEvents(isWellDone,:);
    eventTypes = eventTypes(isWellDone,:);
    nTrials(iFold) = sum(isWellDone);
    
    % Get possibe sequences
    [sequences, chance] = GetAllTrialSequences(regressor_events); %#ok<NASGU> (suppress compiler warning)
    
    % Compare actual data to sequence templates
    iSeq = zeros(1,nTrials(iFold));
    iTrue = zeros(1,nTrials(iFold));
    seqMatch = zeros(nTrials(iFold),size(sequences,1));
    fitMatch = zeros(1,nTrials(iFold));
    isComplete_prediction = false(1,nTrials(iFold));
    for i = 1:nTrials(iFold) % trial
        % Find data for this trial
        iTrialStart = iEvents(i,1)+regressor_range(1)-1;
        iTrialEnd = iEvents(i,end)+regressor_range(2)+1;
        if iTrialEnd>size(EEG.data,2)
            data = [EEG.data(:,iTrialStart:end), nan(size(EEG.data,1),iTrialEnd-size(EEG.data,2))];
        else
            data = EEG.data(:,iTrialStart:iTrialEnd);
        end
        
        % Plot true eeg
        if do_plots
            figure(1); clf; hold on; 
            imagesc(data(chanstoplot,:)); PlotVerticalLines(iEvents(i,:)-iTrialStart,'k--');        
            set(gca,'clim',[-30 30]); colorbar; 
            title(['Fold/Trial ' num2str(iFold) '/' num2str(i) ' TRUE sequence = ' num2str(eventTypes(i,:))]);
            figure(2);
        end
        if do_lines
            figure(3); clf; 
            for k=1:numel(chanstoplot)
                subplot(numel(chanstoplot),1,k); hold on; ylabel(EEG.chanlocs(chanstoplot(k)).labels);
                plot((1:size(data,2))*dt, data(chanstoplot(k),:)/sqrt(var(data(chanstoplot(k),:))),'b'); 
            end
        end
        
        % Get true sequence
        iTrue(i) = find(ismember(sequences,eventTypes(i,:), 'rows'),1);        
        
        % Get EEG-predicted sequence
        for j = 1:size(sequences,1) % sequence
            
            % Build up template (data prediction from this sequence of events) 
            % from sum of offset response functions
            template = zeros(size(data));
            for k = 1:nEventsPerTrial % saccade
                template(:,(iEvents(i,k)-iTrialStart+regressor_range(1)):(iEvents(i,k)-iTrialStart+regressor_range(2))) = ...
                    template(:,(iEvents(i,k)-iTrialStart+regressor_range(1)):(iEvents(i,k)-iTrialStart++regressor_range(2))) + ...
                    responseFns(:,:,sequences(j,k));
            end
            template(template==0) = NaN; % remove time points outside of event influence
                       
            % Get metric for how well this template matches the data
            switch match_method
                case 'AvgRefPearson'
                    matchy = nan(1,size(data,2));
                    for k=1:size(data,2)
                        matchy(k) = corr(data(:,k),template(:,k),'type','Pearson');
                    end
                    seqMatch(i,j) = nanmean(matchy);
                case 'AvgRefNormProduct'
                    avgref_template = template-repmat(mean(template,1),size(template,1),1);
                    seqMatch(i,j) = nanmean(avgref_template(:).*data(:))/sqrt(nanvar(avgref_template(:)));
                case 'Product'
                    seqMatch(i,j) = nanmean(template(:).*data(:));
                case 'NormProduct'
                    seqMatch(i,j) = nanmean(template(:).*data(:))/sqrt(nanvar(template(:)));
                case 'Difference'
                    seqMatch(i,j) = -nanmean(abs(template(:)/sqrt(nanvar(template(:))) - data(:)/sqrt(nanvar(data(:)))));
                otherwise
                    seqMatch(i,j) = corr(data(~isnan(template(:))),template(~isnan(template(:))),'type',match_method);
            end
            % Adjust using the likelihood that this sequence occurs
            if use_posterior
                seqMatch(i,j) = seqMatch(i,j)*chance(j);
            end
            
            % Plot template prediction of true sequence
            if do_lines && j==iTrue(i) 
                for k=1:numel(chanstoplot)
                    subplot(numel(chanstoplot),1,k); hold on;
                    plot((1:size(data,2))*dt, template(chanstoplot(k),:)/sqrt(nanvar(template(chanstoplot(k),:))),'c');
                end
            end
            % Plot template prediction of every sequence
            if do_plots
                subplot(4,4,j); cla; hold on;
                imagesc(template(chanstoplot,:)); PlotVerticalLines(iEvents(i,:)-iTrialStart,'k--');
                set(gca,'clim',[-10 10]); colorbar; title(['sequence = ' num2str(sequences(j,:)) ', matchScore = ' num2str(seqMatch(i,j))]);            
            end
        end
        
        % Determine match of each fit
        [fitMatch(i),iSeq(i)] = max(seqMatch(i,:));
        if pooled_complete
            iCompl = find(strcmp('Compl',regressor_events));
            isComplete = any(sequences==iCompl,2);
            fitMatchComplete = sum(seqMatch(i,isComplete));
            fitMatchIncomplete = sum(seqMatch(i,~isComplete));
            isComplete_prediction(i) = fitMatchComplete>fitMatchIncomplete;
        end
        
        % Label plots of true and EEG-predicted "winner" sequences
        if do_plots
            subplot(4,4,iSeq(i)); text(0,3,' *WINNER*')
            subplot(4,4,iTrue(i)); text(0,1,' *TRUE*')
        end
        % Plot template from EEG-predicted sequence
        if do_lines
            template = zeros(size(data));
            for k = 1:nEventsPerTrial % saccade
                template(:,(iEvents(i,k)-iTrialStart+regressor_range(1)):(iEvents(i,k)-iTrialStart+regressor_range(2))) = ...
                    template(:,(iEvents(i,k)-iTrialStart+regressor_range(1)):(iEvents(i,k)-iTrialStart+regressor_range(2))) + ...
                    responseFns(:,:,sequences(iSeq(i),k));
            end
            template(template==0) = NaN; % remove time points outside of event influence
            % Plot
            for k=1:numel(chanstoplot)
                subplot(numel(chanstoplot),1,k)
                plot((1:size(data,2))*dt, template(chanstoplot(k),:)/sqrt(nanvar(template(chanstoplot(k),:))),'r');
                PlotVerticalLines((iEvents(i,:)-iTrialStart)*dt,'k--');            
            end
            % Annotate plot
            subplot(numel(chanstoplot),1,1); hold on;
            title(sprintf('Fold/Trial %d/%d,  TRUE sequence = %s\n PREDICTED sequence = %s',...
                iFold,i,num2str(eventTypes(i,:)),num2str(sequences(iSeq(i),:))));
            legend('data','template','predicted','events');
            
            % TEMP BARPLOT
            figure(2)
            cla; hold on;           
            isComplete = any(sequences==4,2);
            bar(find(isComplete),seqMatch(i,isComplete),'r');
            bar(find(~isComplete),seqMatch(i,~isComplete),'b');
            plot([iTrue(i) iTrue(i)],[0 seqMatch(i,iTrue(i))],'g-*');
            set(gca,'xtick',1:size(seqMatch,1));
            ylabel(sprintf('seqMatch score (%s)',match_method))
            xlabel('sequence type')
            legend('Complete','Incomplete','True')
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
    prediction{iFold} = sequences(iSeq,:);
    truth{iFold} = eventTypes;
    iCompl = find(strcmp('Compl',regressor_events));
    if pooled_complete
        accuracy(iFold) = mean(isComplete_prediction==any(truth{iFold}==iCompl,2)');
    else
        accuracy(iFold) = mean(any(prediction{iFold}==iCompl,2)==any(truth{iFold}==iCompl,2));
    end    
    seq_accuracy(iFold) = mean(iSeq==iTrue);
    % Display metrics of success
    fprintf('class accuracy: %d/%d = %0.1f%%\n',round(accuracy(iFold)*nTrials(iFold)), nTrials(iFold), accuracy(iFold)*100);
    fprintf('sequence accuracy: %d/%d = %0.1f%%\n',round(seq_accuracy(iFold)*nTrials(iFold)), nTrials(iFold), seq_accuracy(iFold)*100);
end

% Display total accuracy metrics
totalAccuracy = accuracy*nTrials'/sum(nTrials);
totalSeqAccuracy = seq_accuracy*nTrials'/sum(nTrials);
fprintf('TOTAL class accuracy: %d/%d = %0.1f%%\n',round(accuracy*nTrials'), sum(nTrials), totalAccuracy*100);
fprintf('TOTAL sequence accuracy: %d/%d = %0.1f%%\n',round(seq_accuracy*nTrials'), sum(nTrials), totalSeqAccuracy*100);