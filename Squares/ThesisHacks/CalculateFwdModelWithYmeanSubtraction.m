% CalculateFwdModelWithYmeanSubtraction()
% Created 4/19/14 by DJ.

fwdModelNew = zeros(EEG.nbchan,numel(R_train));
for i=1:numel(R_train)
%     if ~isempty(y{i}), continue; end
    fprintf('--- Subject %d/%d...\n',i,numel(R_train));
    % Subtract out other factors
    EEG = R_train(i).EEG;
    EEG = SubtractOutGlmResponses(EEG,R_train(i).responseFns{1},R_train(i).influence{1},R_train(i).regressor_events{1});
    EEG = SubtractOutGlmResponses(EEG,R_train(i).responseFns{2},R_train(i).influence{1},R_train(i).regressor_events{2});  
    
    % Extract epochs   
    EEG0 = pop_epoch( EEG, train_events(1), [-0.5 1], 'newname', train_events{1}, 'epochinfo', 'yes');    
    EEG1 = pop_epoch( EEG, train_events(2), [-0.5 1], 'newname', train_events{2}, 'epochinfo', 'yes');
    
    % Convert twl and two to indices
    t_ms = EEG0.times;
    twl = twl_ms/1000*EEG0.srate; % get window width in samples
    two = round(interp1(t_ms,1:length(t_ms),two_ms)); % get nearest indices
    
    % Assemble data & truth vectors
    data = cat(3,EEG0.data,EEG1.data);
%     truth{i} = [zeros(EEG0.trials,1); ones(EEG1.trials,1)]';
    
    % Run classifier
%     [y{i},w{i},v{i},fwdModel{i},y_level1{i}] = RunHybridHdcaClassifier_LR(data,truth{i},twl,two,'10fold');
%     [y2{i},w2{i},v2{i},fwdModel2{i},y_level12{i}] = RunHybridHdcaClassifier(data,truth{i},twl,two,'10fold');
    meanw = mean(mean(w{i},3),4);
    iTimes = two-1+(1:twl);
    fmData = data(:,iTimes,:);
    fmDataAvg = permute(mean(fmData,2), [3,1,2]);
    yAvg = fmDataAvg*meanw;
    % de-mean y values
    yAvg = yAvg-mean(yAvg);
    % Compute forward model
    fwdModelNew(:,i) = yAvg \ fmDataAvg;


%     figure;
%     PlotHybridHdcaClassifier(fwdModel{i},v{i},EEG0.chanlocs,two_ms + twl_ms/2);
end