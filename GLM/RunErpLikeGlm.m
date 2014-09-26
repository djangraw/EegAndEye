function RunErpLikeGlm(dataset,regressor_events,suffixes,params)

% RunErpLikeGlm(dataset,regressor_events,suffixes,params)
%
% INPUTS:
% -dataset is a string indicating the EEG dataset you want to use as input.
% It should be of format '<prefix>-<subject>-<rest-of-filename>.set'.
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
% Created 3/20/13 by DJ based on RunGlmLikeGui
% Updated 4/3/13 by DJ - fixed SubtractOutGlmResponses inputs
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 3/19/14 by DJ - tResponse, influence, artifact_influence in cells

if nargin<4 || isempty(params)
    params.offset = 0;
    params.influence = [0 500];
    params.artifact_influence = [-500 500];
    params.stddev = 0;
    params.vthresh = 75;
    params.method = 'mvregress';
    params.trial_rej_rules = {'backward','skip','skipped_ends'};
    params.artifact_events = {'Cross','Errant','BlinkStart','BlinkEnd','Button'};
end

if ~isfield(params,'prefix')
    params.prefix = 'sf';
end

% extract subject number
dashes = strfind(dataset,'-');
subject = str2double(dataset((dashes(1)+1):(dashes(2)-1)));

% Clear status
fprintf('%s - Loading input parameters...\n',datestr(now,16));

% Get regressor events
nLevels = length(regressor_events);

% Get filenames
fileprefix = sprintf('%s-%d-ERPresults-',params.prefix,subject);
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

% Turn influence into 2-element vector (tMin, tMax) (NEW 3/20/13)
if numel(influence)==1
    influence = [-abs(influence) abs(influence)];
end
if ~iscell(params.influence)
    params.influence = repmat({params.influence},1,length(regressor_events));
end
if ~iscell(params.artifact_influence)
    params.artifact_influence = repmat({params.artifact_influence},1,length(regressor_events));
end

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


% find events
nEvents = numel(all_events);
nSesssions = numel(y);
[times, codes, weights] = deal(cell(nEvents,nSesssions)); % events by sessions
for i=1:nEvents
    if ~any(strcmp(all_events{i},{EEG.event.type})) % if this event hasn't been added yet...
        [times(i,:),codes(i,:),weights(i,:)] = UseEventRule(y,all_events{i});        
        for j=1:nSesssions
            times{i,j} = times{i,j} + offset; % add offset to times
        end        
    end
end
% append events for each session
[all_times, all_codes] = deal(cell(1,nSesssions));
for j=1:nSesssions
    all_times{j} = cat(1,times{:,j});
    all_codes{j} = cat(1,codes{:,j});    
end

% Add events
EEG = AddEeglabEvents_MultiSession(EEG,y,all_times,all_codes);

% Once all the events are added, add the weights
if ~isfield(EEG.etc,'ureventweights')
    EEG.etc.ureventweights = ones(size(EEG.urevent));
end
for i=1:numel(all_events)
    if ~isempty(weights{i})
        EEG.etc.ureventweights(strcmp(all_events{i},{EEG.urevent(:).type})) = cat(1,weights{i,:});
    end
end

% Get max bounds for influence of events and artifacts
allInfluence = cat(1,influence{:});
maxInfluence = [min(allInfluence(:,1)), max(allInfluence(:,2))]; % max influence interval (at any analysis level)
allArtifact_influence = cat(1,artifact_influence{:});
maxArtifact_influence = [min(allArtifact_influence(:,1)), max(allArtifact_influence(:,2))]; % max influence interval (at any analysis level)

% Epoch data
EEG = pop_epoch(EEG,{'TrialStart-T', 'TrialStart-D'},[-1.5 4.5]);

% Interpolate noisy electrodes and reject noisy trials
EEG = EnforceVoltageThreshold(EEG,vthresh,[regressor_events{:}], artifact_events, maxInfluence, maxArtifact_influence);

% Reject bad trials
EEG = RejectEegData(EEG,y,trial_rej_rules);

% Reject bad events (NEW 3/20/13)
EEG = RedactArtifacts_events(EEG, artifact_events, maxArtifact_influence); % remove any events whose erps would include this artifact.

% Run analysis
for iLevel=1:nLevels % for each level of analysis
    % Update status
    fprintf('%s - Running Analysis Level %d/%d...\n',datestr(now,16),iLevel,nLevels);

    % Get ERP
    erp = zeros(EEG.nbchan,floor(diff(influence{iLevel})/1000*EEG.srate),numel(regressor_events{iLevel}));
    for j=1:numel(regressor_events{iLevel})
        fprintf('level %d, event %d: %s...\n',iLevel,j,regressor_events{iLevel}{j});
        EEGtemp = pop_epoch(EEG,regressor_events{iLevel}(j),influence{iLevel}/1000);
        erp(:,:,j) = mean(EEGtemp.data,3);                
    end        
    
    % Dish out to variables to be saved
    responseFns{iLevel} = erp;
    tResponse{iLevel} = EEGtemp.times;    
    
    if ~isempty(filenames{iLevel})
        % Update status
        fprintf('%s Saving Results as %s...\n',datestr(now,16), filenames{iLevel});        
        % Save results
        save(filenames{iLevel},'EEG','responseFns','tResponse','regressor_events',...
            'filenames','iLevel','artifact_events','dataset','offset','influence',...
            'artifact_influence','stddev','vthresh','method','trial_rej_rules');
    end
    % Subtract out result
    disp('Subtracting out erp...')
%     tResponse_symmetric = -max(abs(tResponse)):(1000/EEG.srate):max(abs(tResponse));
%     erp_symmetric = zeros(size(erp,1),length(tResponse_symmetric),size(erp,3));
%     erp_symmetric(:,tResponse_symmetric>=tResponse(1) & tResponse_symmetric<=tResponse(end),:) = erp;            
%     EEG = SubtractOutGlmResponses(EEG,erp_symmetric,regressor_events{iLevel});
    EEG = SubtractOutGlmResponses(EEG,erp,[tResponse(1) tResponse(end)],regressor_events{iLevel});
end

% Update status
fprintf('%s - Analysis complete!\n',datestr(now,16));