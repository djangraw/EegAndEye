function RunAddOnLevelOfGlm(results,baseLevel,new_regressor_events,new_suffixes)

% RunAddOnLevelOfGlm(results,baseLevel,regressor_events,suffixes)
%
% INPUTS:
% -results is a string indicating the results dataset you want to start
% with, or a struct of that dataset already loaded.
% -baseLevel is a scalar indicating the level of results that you want to
% start after.
% -regressor_events is an N-element cell array of cell arrays of strings
% indicating what events should be used in your N additional levels of GLM.
% e.g., {{'LSaccade' 'RSaccade'} {'SqNum1' 'SqNum2' 'SqNum3'}}.
% -suffixes is an N-element cell array of strings indicating what the
% output of each GLM level should be saved as.
% e.g., {'LRSaccade-v1pt0.mat' 'LRSaccade-SqNum-v1pt0.mat'}.
%
% Other parameters will be taken from the results dataset.
%
% Created 1/24/13 by DJ- adapted from RunGlmGui.m "Run" button
% Updated 4/4/13 by DJ - fixed artifact_influence bug
% Updated 4/15/13 by DJ - added SqFix compatibility
% Updated 4/30/13 by DJ - made responseFns a cell array
% Updated 3/19/14 by DJ - tResponse in cells

% Unpack variables from results struct
if ischar(results)
    fprintf('loading %s...\n',results);
    results = load(results);
end
results = UpdateGlmResultsFormat(results);
if isstruct(results)
    fprintf('unpacking results struct...\n');
    UnpackStruct(results);
end
if ~exist('dataset','var'), dataset = []; end % trick MATLAB compiler, b/c dataset is a matlab fn.
if ~exist('artifact_influence','var'), artifact_influence = influence; end
    
% Get subject number
dashes = strfind(dataset,'-');
dataprefix = dataset(1:(dashes(1)-1));
subject = str2double(dataset((dashes(1)+1):(dashes(2)-1)));

% Clear status
fprintf('%s - Loading input parameters...\n',datestr(now,16));

% Get regressor events
% baseLevel = iLevel;
regressor_events = [regressor_events(1:baseLevel); new_regressor_events];
nLevels = length(regressor_events);

% Get filenames
prefix = sprintf('%s-%d-GLMresults-',dataprefix,subject);
new_filenames = cell(size(new_suffixes));
for i=1:numel(new_suffixes)
    if isempty(new_suffixes{i})
        new_filenames{i} = '';
    else
        new_filenames{i} = [prefix, new_suffixes{i}];
    end
end
filenames = [filenames(1:baseLevel); new_filenames];

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
y = loadBehaviorData(subject,[],dataprefix);
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
        EEG.etc.ureventweights(strcmp(all_events{i},{EEG.urevent(:).type})) = cat(1,weights{i,:})';
    end
end


% Epoch data
EEG = pop_epoch(EEG,{'TrialStart-T', 'TrialStart-D'},[-1.5 4.5]);

% Interpolate noisy electrodes and reject noisy trials
EEG = EnforceVoltageThreshold(EEG,vthresh,[regressor_events{:}], artifact_events, influence, artifact_influence);

% Reject bad trials
EEG = RejectEegData(EEG,y,trial_rej_rules);

% Subtract out early level response functions
EEGnew = EEG; % create copy to be modified during run (but only original will be saved)
for iLevel = 1:baseLevel
    if isempty(responseFns{iLevel}) % For older datasets
        % Update status
        fprintf('%s - Loading Analysis Level %d/%d: %s...\n',datestr(now,16),iLevel,nLevels,filenames{iLevel});
        R = load(filenames{iLevel});
        R = UpdateGlmResultsFormat(R);
        fprintf('%s - Subtracting out Analysis Level %d/%d...\n',datestr(now,16),iLevel,nLevels);
        responseFns{iLevel} = R.responseFns{iLevel};
        EEGnew = SubtractOutGlmResponses(EEGnew,R.responseFns{iLevel},R.influence,regressor_events{iLevel});    
    else % 4/30/13 version
        % Update status        
        fprintf('%s - Subtracting out Analysis Level %d/%d...\n',datestr(now,16),iLevel,nLevels);        
        EEGnew = SubtractOutGlmResponses(EEGnew,responseFns{iLevel},influence,regressor_events{iLevel});    
    end
end
    


% Run analysis
for iLevel = (baseLevel+1):nLevels % for each level of analysis
    % Update status
    fprintf('%s - Running Analysis Level %d/%d...\n',datestr(now,16),iLevel,nLevels);

    [responseFns{i}, tResponse{iLevel}] = RunEegGlm(EEGnew,regressor_events{iLevel},influence,artifact_events,artifact_influence,method,stddev);
%     responseFns = [responseFns; {newRF}];
%     responseFns = []; tResponse  = [];
    
    if ~isempty(filenames{iLevel})
        % Update status
        fprintf('%s Saving Results as %s...\n',datestr(now,16), filenames{iLevel});
        % Save results
        save(filenames{iLevel},'EEG','responseFns','tResponse','regressor_events',...
            'filenames','iLevel','artifact_events','artifact_influence','dataset','offset','influence','stddev',...
            'vthresh','method','trial_rej_rules');
    end
    EEGnew = SubtractOutGlmResponses(EEGnew,responseFns{iLevel},influence,regressor_events{iLevel});
end

% Update status
fprintf('%s - Analysis complete!\n',datestr(now,16));