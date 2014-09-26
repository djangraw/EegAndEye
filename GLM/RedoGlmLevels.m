function RedoGlmLevels(R,newParams)

% RedoGlmLevels(R,newParams)
%
% Created 3/19/14 by DJ.
% Updated 8/14/14 by DJ - added lambda input to RunEegGlm.

% Update results struct
R = UpdateGlmResultsFormat(R);

UnpackStruct(R);
UnpackStruct(newParams); % overwrite any params you want to change
nLevels = length(regressor_events);
% Fool compiler by declaring dummy variables
if ~exist('dataset','var'), dataset = ''; end
if ~exist('trial_rej_rules','var'), trial_rej_rules = ''; end
if ~exist('lambda','var'),lambda = NaN; end
if ~exist('demeanX','var'),demeanX = false; end

% Check for match
oldevents = [R.artifact_events, R.regressor_events{:}];
newevents = [artifact_events, regressor_events{:}];
if ~isempty(setdiff(newevents,oldevents)) % if any newevents were not in oldevents
    error('Some events not in dataset!');
end
if ~strcmp(dataset,R.dataset) || ~isequal(trial_rej_rules,R.trial_rej_rules)
    error('Datasets do not match!');
end

% Subtract out early level response functions
EEGnew = EEG; % create copy to be modified during run (but only original will be saved)
matchSoFar = true; % have all the levels matched until now?
for iLevel=1:nLevels
    % if it's a match, subtract out existing results
    if matchSoFar && isequal(regressor_events{iLevel},R.regressor_events{iLevel}) && isequal(influence{iLevel}, R.influence{iLevel}) && isequalwithequalnans(lambda{iLevel}, R.lambda{iLevel})
        % Update status  
        fprintf('Level %d matches previous analysis.\n',iLevel);
        fprintf('%s - Subtracting out Analysis Level %d/%d...\n',datestr(now,16),iLevel,nLevels);        
        EEGnew = SubtractOutGlmResponses(EEGnew,R.responseFns{iLevel},R.influence{iLevel},R.regressor_events{iLevel});            
    else
        matchSoFar = false;
        responseFns = responseFns(1:iLevel-1);
        tResponse = tResponse(1:iLevel-1);
        % Update status
        fprintf('Level %d does not match previous analysis.\n',iLevel);
        fprintf('%s - Running Analysis Level %d/%d...\n',datestr(now,16),iLevel,nLevels);

        [newRF, newTR,~,lambda{iLevel}] = RunEegGlm(EEGnew,regressor_events{iLevel},influence{iLevel},artifact_events,artifact_influence{iLevel},method,stddev,lambda{iLevel},demeanX);
        responseFns{iLevel} = newRF;
        tResponse{iLevel} = newTR;
    %     responseFns = []; tResponse  = [];

        if ~isempty(filenames{iLevel})
            % Update status
            fprintf('%s Saving Results as %s...\n',datestr(now,16), filenames{iLevel});
            % Save results
            save(filenames{iLevel},'EEG','responseFns','tResponse','regressor_events',...
                'filenames','iLevel','artifact_events','artifact_influence','dataset','offset','influence','stddev',...
                'vthresh','method','trial_rej_rules','lambda','demeanX');
        end
        EEGnew = SubtractOutGlmResponses(EEGnew,responseFns{iLevel},influence{iLevel},regressor_events{iLevel});
    end
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