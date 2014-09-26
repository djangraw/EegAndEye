function RunXfoldEegGlm(EEG, nuisance_events, regressor_events, artifact_events, extent_ms, artifact_extent_ms, reg_method, nFolds)

% RunXfoldEegGlm(EEG, nuisance_events, regressor_events, artifact_events, 
%   extent_ms, reg_method, nFolds)
%
% INPUTS:
% -EEG is an eeglab dataset that has been run through SetUpGlm
% -nuisance_events, regressor_events, and artifact_events are the event
% name outputs of SetUpGlm
% -extent_ms is a scalar indicating the extent of each event's influence,
% in ms
% -reg_method is a string indicating the regression method you want to use
% in the GLM ('mvregress' is typical')
% -nFolds is a scalar indicating the number of folds you want to use.
%
% OUTPUTS: saves these files in the current directory
% -sq-<subject>-GLMresults-Nuisance-fold<1-nFolds>
% -sq-<subject>-GLMresults-Nuisance-Training-fold<1-nFolds>
% -sq-<subject>-GLMresults-Nuisance-Testing-fold<1-nFolds>
%
% Created 2/1/12 by DJ.
% Updated 3/19/12 by DJ - suppress code warnings
% Updated 3/22/13 by DJ - added artifact_extent_ms input
% Updated 4/30/13 by DJ - responseFns in cells

% Set up
hyphens = strfind(EEG.setname,'-');
subject = str2double(EEG.setname( (hyphens(1)+1):(hyphens(2)-1) ));
foldBounds = GetGlmFoldBounds(EEG,regressor_events,artifact_events,extent_ms,nFolds,artifact_extent_ms);

% Run folds
for i=1:nFolds
    % Separate out training and testing sets
    fprintf('---Fold %d of %d...\n',i,nFolds);
    EEGtrain_init = eeg_eegrej(EEG, foldBounds([i i+1]));
    EEGtest_init = eeg_eegrej(EEG, [foldBounds([1 i]); foldBounds([i+1 end])]);
    
    % Subtract out nuisance events
    if ~isempty(nuisance_events)
        % Run nuisance GLM
        disp('Running Nuisance GLM...');
        [newRF, tResponse, NewEEG] = RunEegGlm(EEGtrain_init,...
            nuisance_events, extent_ms, artifact_events,...
            artifact_extent_ms,reg_method);
        responseFns{1} = newRF;
        
        disp('Saving...')
        save(sprintf('sq-%d-GLMresults-Nuisance-fold%d',subject,i),...
            'responseFns','tResponse','NewEEG','nuisance_events','artifact_events','extent_ms','artifact_extent_ms');
        
        % Subtract nuisance regressors out of both training and testing data
        disp('Subtracting...')
        EEGtrain = SubtractOutGlmResponses(EEGtrain_init,responseFns{1},artifact_extent_ms,nuisance_events);
        EEGtest = SubtractOutGlmResponses(EEGtest_init,responseFns{1},artifact_extent_ms,nuisance_events);
    else
        % Skip nuisance regressor
        EEGtrain = EEGtrain_init;
        EEGtest = EEGtest_init;
    end
    
    % Run Training GLM
    disp('Running Training GLM...')
    [newRF, tResponse, NewEEG] = RunEegGlm(EEGtrain,regressor_events,extent_ms,artifact_events,artifact_extent_ms,reg_method); %#ok<*NASGU,*ASGLU> (suppresses code warnings)
    responseFns{2} = newRF;
    
    % Save results
    disp('Saving...')
    save(sprintf('sq-%d-GLMresults-Nuisance-Training-fold%d',subject,i),...
        'responseFns','tResponse','NewEEG','nuisance_events','regressor_events','artifact_events','extent_ms','artifact_extent_ms');
    
    % Prep test data
    results = load(sprintf('sq-%d-GLMresults-Nuisance-Training-fold%d',subject,i));
    [h, y, iEvents, truth] = PrepareGlmDataForClassification(results,EEGtest);
    save(sprintf('sq-%d-EEGtest-Nuisance-Testing-fold%d',subject,i),'EEGtest','h','y','iEvents','truth','regressor_events');
    disp('Done!')
end