function [NewEEG,S] = SubtractOutGlmResponses(EEG,responseFns,extent_ms,event_types,stddev,normalizeX,demeanX)

% Subtract GLM responses from EEG data each time specified events occur.
%
% [NewEEG,S] = SubtractOutGlmResponses(EEG,responseFns,extent_ms,event_types,stddev,normalizeX,demeanX)
%
% Note that we choose not to address artifact events here.
%
% INPUTS:
% -EEG is an eeglab data struct with the specified event_types added.
% -responseFns is a DxTxN struct, where D is the number of channels in the
%  EEG, T is the number of time points, and N is the number of event types
%  you want to subtract out.
% -extent_ms is a scalar or 2-element vector indicating the time relative to the
% event that responseFns starts and ends.
% -event_types is a 1xN cell array of strings indicating the event labels
%  corresponding to each response function.
% -stddev is a scalar indicating the width of the gaussian you want to
% convolve with the regressor matrix (default = 0).
%
% OUTPUTS:
% -NewEEG is an eeglab data struct identical to EEG except the data has had
%  the requested responses subtracted out.
% -S is the regressor model matrix.
%
% Created 1/12/12 by DJ.
% Updated 1/31/12 by DJ - 'event not found' check
% Updated 8/8/12 by DJ - added support for epoched datasets
% Updated 8/9/12 by DJ - added support for event_weights, stddev
% Updated 3/22/13 by DJ - added input for extent_ms, GetGlmRegressors_v2p0.
% Updated 7/31/14 by DJ - added S output.

disp('Setting up...')
if nargin<5 || isempty(stddev)
    stddev = 0;
end
if nargin<6 || isempty(normalizeX)
    normalizeX = false;
end
if nargin<7 || isempty(demeanX)
    demeanX = false;
end

% Find constants
dt = 1000/EEG.srate;
regressor_range = round(extent_ms/dt); % how many samples should each response function extend?
% Nt = (size(responseFns,2)-1)/2;
t = (1:EEG.pnts*EEG.trials)*dt; % time vector for EEG

% Find event times
Nr = numel(event_types);
event_times = cell(1,Nr);
event_weights = cell(1,Nr);
for i=1:Nr
    isGoodEvent = strcmp(event_types{i},{EEG.event(:).type}); % Don't worry about trial rejection
    event_times{i} = [EEG.event(isGoodEvent).latency]*dt;
    if isfield(EEG.etc,'ureventweights') % v2pt1 results and after
        event_weights{i} = EEG.etc.ureventweights([EEG.event(isGoodEvent).urevent]);
    elseif isfield(EEG.etc,'eventweights') % v2pt0 results
        event_weights{i} = EEG.etc.eventweights(isGoodEvent);
    else % v1ptX results
        event_weights{i} = ones(size(event_times{i}));
    end
    % Warning if event is not found
    if isempty(event_times{i})
        warning('DJ:SubtractOutGlmResponses:EventNotFound','Event %s not found!\n',event_types{i});
    end
end

disp('Getting regressors...');
% Get regressors
% [~,S] = GetGlmRegressors(t,event_times,[],regressor_range,event_weights,stddev); % Don't worry about artifact events
[~,S] = GetGlmRegressors_v2p0(t,event_times,[],regressor_range,0,event_weights,stddev); % Don't worry about artifact events

% DE-MEAN AND/OR NORMALIZE S MATRIX
if normalizeX
    S = full(S);
    isInPlay = any(S~=0,2);
    X = full(S(isInPlay,:));
    if demeanX
        fprintf('De-meaning and ');
        mnX = mean(X,1);
    else
        fprintf('NOT de-meaning, but ');
        mnX = zeros(1,size(X,2));
    end
    fprintf('Standardizing design matrix...\n')
    p = size(X,2);
    rssX = nan(size(mnX));
    for i=1:p  
        foo = (X(:,i) - mnX(i));
        rssX(i) = sqrt(sum(foo.^2)); % root sum square
        X(:,i) = foo / rssX(i);
    end 
    S(isInPlay,:) = X;
end




disp('Subtracting out response...');
% Multiply response functions by regressors
foo = permute(responseFns,[2 3 1]);
H = reshape(foo,[size(foo,1)*size(foo,2),size(foo,3)]);
reconstructed = (S*H)';

NewEEG = EEG;
NewData = reshape(EEG.data,EEG.nbchan,EEG.pnts*EEG.trials) - reconstructed; % residuals
NewEEG.data = reshape(NewData,size(EEG.data));

disp('Succcess!')