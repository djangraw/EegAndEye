function [X,Y,Xmean,Xrss] = GetGlmMatrices(EEG,event_types,extent_ms,artifact_types,artifact_extent_ms,demeanX,normalizeX)

% [X,Y] = GetGlmMatrices(EEG,event_types,extent_ms,artifact_types,artifact_extent_ms,demeanX,normalizeX)
%
% INPUTS:
%
% OUTPUTS:
% -X is an nxp regressor matrix.
% -Y is an nxD matrix of observed data.
%
% Created 8/27/14 by DJ.


% Find constants
dt = 1000/EEG.srate;
regressor_range = round(extent_ms/dt); % how many samples should each response function extend?
artifact_range = round(artifact_extent_ms/dt); % how many samples should each artifact affect?
t = (1:EEG.pnts*EEG.trials)*dt; % time vector for EEG

% Find event times
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
artifact_times = [EEG.event(ismember({EEG.event(:).type}, artifact_types)).latency]*dt;

disp('Getting regressors...');
% Get regressors
[~,S] = GetGlmRegressors_v2p0(t,event_times,artifact_times,regressor_range,artifact_range,event_weights);

% Crop data
isRelevantTime = sum(S,2)>0;
X = full(S(isRelevantTime,:));
% tcrop = t(isRelevantTime);
X = double(X);
if demeanX
    disp('De-meaning X...')
    Xmean = mean(X,1);
    for i=1:size(X,2)        
        X(:,i) = X(:,i)-Xmean(i);
    end
else
    Xmean = zeros(1,size(X,2));
end
if normalizeX
    disp('Normalizing X...')
    Xrss = ones(1,size(X,2));
    for i=1:size(X,2)
        Xrss(i) = sqrt(X(:,i)'*X(:,i));
        X(:,i) = X(:,i)/Xrss(i);
    end
else
    Xrss = ones(1,size(X,2));
end

if nargout>1
    % Reshape data, if necessary
    reshapeddata = reshape(EEG.data,size(EEG.data,1),size(EEG.data,2)*size(EEG.data,3))';
    % Crop reshaped data
    Y = reshapeddata(isRelevantTime,:);
end
