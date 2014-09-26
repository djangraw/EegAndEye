function [XXinvX,Y,X] = GetGlmInverseMatrix(R,lambda)

% Created 8/8/14 by DJ.

UnpackStruct(R);
method = 'ridge'; % overwrite method leastsquares

% Reshape data, if necessary
reshapeddata = reshape(EEG.data,size(EEG.data,1),size(EEG.data,2)*size(EEG.data,3));

% Find constants
dt = 1000/EEG.srate;
regressor_range = round(influence{iLevel}/dt); % how many samples should each response function extend?
artifact_range = round(artifact_influence{iLevel}/dt); % how many samples should each artifact affect?
t = (1:size(reshapeddata,2))*dt; % time vector for EEG

% Find event times
event_types = regressor_events{iLevel};
Nr = numel(event_types);
event_times = cell(1,Nr);
event_weights = cell(1,Nr);
if ~isfield(EEG.etc,'rejectepoch') % Add rejectepoch field if it's not there already
    EEG.etc.rejectepoch = zeros(EEG.trials,1);
end
for i=1:Nr
    isGoodEvent = strcmp(event_types{i},{EEG.event(:).type}) & ~EEG.etc.rejectepoch([EEG.event(:).epoch])';
    event_times{i} = [EEG.event(isGoodEvent).latency]*dt;
    event_weights{i} = EEG.etc.ureventweights([EEG.event(isGoodEvent).urevent]);
end

% Find blink events
artifact_times = [EEG.event(ismember({EEG.event(:).type}, artifact_events)).latency]*dt;

disp('Getting regressors...');
% Get regressors
% [~,S] = GetGlmRegressors(t,event_times,artifact_times,Nt,event_weights,stddev);
[~,S] = GetGlmRegressors_v2p0(t,event_times,artifact_times,regressor_range,artifact_range,event_weights,stddev);

% Crop data
isRelevantTime = sum(S,2)>0;
Scrop = S(isRelevantTime,:);
Ecrop = reshapeddata(:,isRelevantTime);
% tcrop = t(isRelevantTime);

% prepare outputs
X = full(double(Scrop));
Y = Ecrop';

% Run real regression
disp('Calculating Inverse...');
XXinvX = (X'*X + lambda*eye(size(X,2)))^(-1)*X';
disp('Done!')