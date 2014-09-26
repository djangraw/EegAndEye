function RunGlmLikeGui(dataset,regressor_events,suffixes,params)

% RunGlmLikeGui(dataset,regressor_events,suffixes,params)
%
% INPUTS:
% -dataset is a string indicating the EEG dataset you want to use as input.
% It should be of format 'sq-<subject>-<rest-of-filename>.set'.
% -regressor_events is an N-element cell array of cell arrays of strings
% indicating what events should be used in your N-level GLM.
% e.g., {{'LSaccade' 'RSaccade'} {'SqNum1' 'SqNum2' 'SqNum3'}}.
% -suffixes is an N-element cell array of strings indicating what the
% output of each GLM level should be saved as.
% e.g., {'LRSaccade-v1pt0.mat' 'LRSaccade-SqNum-v1pt0.mat'}.
% -params is a struct with fields:
% offset,influence,stddev,vthresh,method,trial_rej_rules,artifact_events.
% See RunGuiGlm for details.
%
% Created 1/24/13 by DJ- adapted from RunGlmGui.m "Run" button
% Updated 3/6/13 by DJ - added prefix parameter
% Updated 4/23/13 by DJ - messed with event weights to fix a sf-10 problem!
% Updated 4/26/13 by DJ - replaced RejectEegData w/ MatchFixationDurations
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 5/10/13 by DJ - delete lower-level files
% Updated 3/19/14 by DJ - tResponse, influence, and artifact_influence in
%  cells. Calls MatchFixationDurations and EnforceVoltageThreshold using 
%  largest influence bounds.
% Updated 8/14/14 by DJ - added lambda input to RunEegGlm
% Updated 8/26/14 by DJ - added demeanX parameter.

if nargin<4 || isempty(params)
    params.offset = 0;
    params.influence = [-500 500];
    params.artifact_influence = [-500 500];
    params.stddev = 0;
    params.vthresh = 75;
    params.method = 'mvregress';
    params.trial_rej_rules = {'backward','skip','skipped_ends'};
    params.artifact_events = {'Cross','Errant','BlinkStart','BlinkEnd','Button'};
    params.contrast_events = {};
end

if ~isfield(params,'prefix')
    params.prefix = 'sq';
end
if ~iscell(params.influence)
    params.influence = repmat({params.influence},1,length(regressor_events));
end
if ~iscell(params.artifact_influence)
    params.artifact_influence = repmat({params.artifact_influence},1,length(regressor_events));
end
if ~isfield(params,'lambda')
    params.lambda = repmat({NaN},1,length(regressor_events));
end
if ~isfield(params,'demeanX')
    params.demeanX = false;
end

% extract subject number
dashes = strfind(dataset,'-');
subject = str2double(dataset((dashes(1)+1):(dashes(2)-1)));

% Clear status
fprintf('%s - Loading input parameters...\n',datestr(now,16));

% Get regressor events
nLevels = length(regressor_events);

% Get filenames
fileprefix = sprintf('%s-%d-GLMresults-',params.prefix,subject);
filenames = cell(size(suffixes));
for i=1:numel(suffixes)
    if isempty(suffixes{i})
        filenames{i} = '';
    else
        filenames{i} = [fileprefix, suffixes{i}];
    end
end

% Get settings
UnpackStruct(params);

% Update status
fprintf('%s - Preparing dataset %s...\n',datestr(now,16),dataset);

% Load data
eeglab nogui
if strcmp(dataset(end-3:end),'.set')
    EEG = pop_loadset(dataset);
else
    results = load(dataset);
%     EEG = SubtractOutGlmResponses(results.EEG,results.responseFns,results.regressor_events{results.iLevel});
    EEG = results.EEG;
end

% Add specified events as in SetUpGLM
% load behavioral data
y = loadBehaviorData(subject,[],prefix);
% add events
all_events = [artifact_events regressor_events{:}];
% for iLevel=1:nLevels
%     all_events = [all_events, regressor_events{iLevel}];
% end

% find events
nEvents = numel(all_events);
nSessions = numel(y);
[times, codes, weights] = deal(cell(nEvents,nSessions)); % events by sessions
for i=1:nEvents
    if ~any(strcmp(all_events{i},{EEG.event.type})) % if this event hasn't been added yet...
        [times(i,:),codes(i,:),weights(i,:)] = UseEventRule(y,all_events{i});        
        for j=1:nSessions
            times{i,j} = times{i,j} + offset; % add offset to times
        end        
    end
end 

% append events for each session
[all_times, all_codes, all_weights] = deal(cell(1,nSessions));
for j=1:nSessions
    all_times{j} = cat(1,times{:,j});
    all_codes{j} = cat(1,codes{:,j});
    all_weights{j} = cat(1,weights{:,j}); % added 4/23/13
end

% Get max bounds for influence of events and artifacts
allInfluence = cat(1,influence{:});
maxInfluence = [min(allInfluence(:,1)), max(allInfluence(:,2))]; % max influence interval (at any analysis level)
allArtifact_influence = cat(1,artifact_influence{:});
maxArtifact_influence = [min(allArtifact_influence(:,1)), max(allArtifact_influence(:,2))]; % max influence interval (at any analysis level)


% NEW 4/26/13: Use MatchFixationDurations to accept/reject data
if ~isempty(contrast_events)    
    tLong = diff(maxInfluence); 
    nTrialsPerSession = numel(y(1).trial.start_time);
    [isKeeper, class] = MatchFixationDurations(y,contrast_events,tLong,trial_rej_rules);
    for k=1:length(regressor_events)
        [~, class] = MatchFixationDurations(y,regressor_events{k},[],[],'none'); % shortcut to get classes of corresponding events
        for j=1:nSessions
            isK_this_sess = isKeeper((1:nTrialsPerSession)+(j-1)*nTrialsPerSession,:);
            class_this_sess = class((1:nTrialsPerSession)+(j-1)*nTrialsPerSession,:);
            for i=1:nEvents
                iThisCode = find(strcmp(all_events{i},all_codes{j})); % indices within all_...{j}
                is_this_class = strcmp(all_events{i},class_this_sess); % indices within ..._this_sess
                isout = ~isK_this_sess(is_this_class); % all_... indices to be removed
                all_times{j}(iThisCode(isout)) = [];
                all_codes{j}(iThisCode(isout)) = [];
                all_weights{j}(iThisCode(isout)) = [];
            end
        end
    end
end


% Add events
[EEG,~,~,iRemoved] = AddEeglabEvents_MultiSession(EEG,y,all_times,all_codes);
% Remove these events from the lists
for j=1:nSessions
    all_weights{j}(iRemoved{j}) = [];
    all_codes{j}(iRemoved{j}) = [];
    all_times{j}(iRemoved{j}) = [];
end

% Once all the events are added, add the weights
if ~isfield(EEG.etc,'ureventweights')
    EEG.etc.ureventweights = ones(size(EEG.urevent));
end
% EEG.etc.ureventweights(ismember({EEG.urevent(:).type},all_events)) = 1;%cat(1,all_weights{:}); % FIX THIS!!!!
% Replace with true weights
all_weights_vec = cat(1,all_weights{:});
all_codes_vec = cat(1,all_codes{:});
for i=1:numel(all_events)
    these_weights = all_weights_vec(strcmp(all_events{i},all_codes_vec));
    EEG.etc.ureventweights(strcmp(all_events{i},{EEG.urevent(:).type})) = these_weights;    
end


% Epoch data
EEG = pop_epoch(EEG,{'TrialStart-T', 'TrialStart-D'},[-1.5 4.5]);

% Interpolate noisy electrodes and reject noisy trials
EEG = EnforceVoltageThreshold(EEG,vthresh,[regressor_events{:}], artifact_events, maxInfluence, maxArtifact_influence);

% Reject bad trials 
if isempty(contrast_events) % if MatchFixationDuration hasn't done it...
    EEG = RejectEegData(EEG,y,trial_rej_rules);
end

EEGnew = EEG; % create copy to be modified during run (but only original will be saved)
% Run analysis
for iLevel=1:nLevels % for each level of analysis
    % Update status
    fprintf('%s - Running Analysis Level %d/%d...\n',datestr(now,16),iLevel,nLevels);

    [newRF, newTR, ~, lambda{iLevel}] = RunEegGlm(EEGnew,regressor_events{iLevel},influence{iLevel},artifact_events,artifact_influence{iLevel},method,stddev,lambda{iLevel},demeanX);
    responseFns{iLevel} = newRF;
    tResponse{iLevel} = newTR;
%     responseFns = []; tResponse  = [];
    
    if ~isempty(filenames{iLevel})
        % Update status
        fprintf('%s Saving Results as %s...\n',datestr(now,16), filenames{iLevel});        
        % Save results
        save(filenames{iLevel},'EEG','responseFns','tResponse','regressor_events',...
            'filenames','iLevel','artifact_events','dataset','offset','influence',...
            'artifact_influence','stddev','vthresh','method','trial_rej_rules','lambda','demeanX');
    end
    EEGnew = SubtractOutGlmResponses(EEGnew,responseFns{iLevel},influence{iLevel},regressor_events{iLevel},stddev);
end

% Clean up by deleting lower-level files
disp('Cleaning up...')
for iLevel = 1:nLevels-1
    if ~strcmp(filenames{iLevel},filenames{nLevels})
        try
            fprintf('Deleting file %s...\n',filenames{iLevel});
            delete(filenames{iLevel});
        catch
            warning('Couldn''t delete file %s!',filenames{iLevel});
        end
    end
end

% Update status
fprintf('%s - Analysis complete!\n',datestr(now,16));