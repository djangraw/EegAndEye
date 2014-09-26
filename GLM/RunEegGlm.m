function [responseFns, tResponse, NewEEG, lambda] = RunEegGlm(EEG, event_types, extent_ms, artifact_types, artifact_extent_ms, reg_method, stddev,lambda,demeanX)

% Estimates the EEG responses to various events using a general linear
% model (GLM).
%
% [responseFns, tResponse, NewEEG] = RunEegGlm(EEG, event_types, extent_ms,
% artifact_types, artifact_extent_ms, reg_method, stddev)
%
% INPUTS:
% -EEG is an unepoched EEGLAB dataset (size(EEG.data,3) = 1).
% -event_types is a cell array of strings indicating the names of events
%  you want to make into regressors.
% -extent_ms is a scalar or 2-element vector indicating the number of ms 
%  you want to include in each regressor's extent.
% -artifact_types is a cell array of strings indicating the names of events
%  that cause EEG artifacts.  All regressors within extent_ms of these
%  artifact events will be zeroed out.
% -artifact_extent_ms is a scalar or 2-element vector indicating the number of ms 
%  you want to zero out relative to each artifact event.
% -reg_method is a string indicating the method you want to use for
%  multivariate regression.  The default is 'mvregress.'  See SimpleGlm.m
%  for the other options.
% -stddev is a scalar indicating the standard deviation (in samples) of 
%  the gaussian you wish to convolve with the regressor matrix.  Default is
%  0, indicating a delta function at the time of each event.
% -lambda is an optional scalar used as the penalization term in ridge
%  regression (reg_method='ridge'). [default = 1e-6].
% 
% OUTPUTS:
% -responseFns is a DxTxNr matrix, where D is the number of channels in the
%  EEG dataset, T is the number of timepoints included in the regressor, 
%  and Nr is the number of regressors.  responseFns(i,j,k) is the estimated
%  response of electrode i to regressor k at sample j.
% -tResponse is a T-element vector showing the times of the responses
%  relative to the regressors.  EXAMPLE: plot(tResponse,responseFns(i,:,k) 
%  will plot the response of electrode i to regressor k.
% -NewEEG is an EEGLAB data file whose data is that reconstructed from the
%  regressors.
%
% Created 12/22/11 by DJ.
% Updated 12/23/11 by DJ.
% Updated 12/27/11 by DJ - parfor!
% Updated 1/10/12 by DJ - useRidge option
% Updated 1/11/12 by DJ - changed useRidge to reg_method, added try/catch
% Updated 6/18/12 by DJ - added event_weights() variable to handle LogOdds
% events.
% Updated 7/16/12 by DJ - added stddev input
% Updated 8/7/12 by DJ - added reshapeddata for epoched-data support (NOTE:
% regressor events must be at least extent_ms*2 ms from the epoch limits!)
% Updated 8/8/12 by DJ - added EEG.etc.rejectepoch support
% Updated 3/22/13 by DJ - use asymmetric influence, add artifact_extent_ms
% input, use GetGlmRegressors_v2p0
% Updated 4/22/13 by DJ - fixed rejectepoch default
% Updated 8/3/14 by DJ - added leastsquares option
% Updated 8/8/14 by DJ - added ridge option, lambda input
% Updated 8/26/14 by DJ - added demeanX option.

if nargin<6 || isempty(reg_method)
    reg_method = '';
end
if nargin<7 || isempty(stddev)
    stddev = 0;
end
if nargin<8 || isempty(lambda)
    lambda = NaN;
end
if nargin<9 || isempty(demeanX)
    demeanX = false;
end

% Handle 1-element range
if numel(extent_ms)==1
    extent_ms = [-extent_ms extent_ms];
end


disp('Setting up...');
% Set options
plotResults = 0;
chanToPlot = 'AF7';

% Reshape data, if necessary
reshapeddata = reshape(EEG.data,size(EEG.data,1),size(EEG.data,2)*size(EEG.data,3));

% Find constants
dt = 1000/EEG.srate;
regressor_range = round(extent_ms/dt); % how many samples should each response function extend?
artifact_range = round(artifact_extent_ms/dt); % how many samples should each artifact affect?
t = (1:size(reshapeddata,2))*dt; % time vector for EEG

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
%     wts =  GetEventWeights(str2double(EEG.subject),event_types{i}); % get weights (e.g., for LogOdds events.)  For unweighted events, will be blank.
%     event_weights{i} = cat(2,wts{:}); % concatenate across all sessions
end

% Find blink events
artifact_times = [EEG.event(ismember({EEG.event(:).type}, artifact_types)).latency]*dt;

disp('Getting regressors...');
% Get regressors
% [~,S] = GetGlmRegressors(t,event_times,artifact_times,Nt,event_weights,stddev);
[~,S] = GetGlmRegressors_v2p0(t,event_times,artifact_times,regressor_range,artifact_range,event_weights,stddev);

% Crop data
isRelevantTime = sum(S,2)>0;
Scrop = S(isRelevantTime,:);
Ecrop = reshapeddata(:,isRelevantTime);
tcrop = t(isRelevantTime);

disp('Reconstructing data...');
% Reconstruct data
NewData = cell(1,EEG.nbchan);
% NewData = zeros(size(EEG.data));
responseFns = nan([Nr, diff(regressor_range)+1, EEG.nbchan]);
chanLabels = {EEG.chanlocs(:).labels}; % to speed up computation

if strcmp(reg_method,'leastsquares')
    responseFns = SimpleGLM(Ecrop', tcrop, double(Scrop), regressor_range, reg_method);
elseif strcmp(reg_method,'ridge')
    [responseFns,~,~,~,lambda] = SimpleGLM(Ecrop', tcrop, double(Scrop), regressor_range, reg_method,lambda,demeanX);
else
    RF = cell(1,EEG.nbchan);
    parfor i=1:EEG.nbchan
        tic;
        fprintf('Electrode %s...',chanLabels{i})
    %     [responseFns(:,:,i),~,tmp] = SimpleGLM(Ecrop(i,:)', tcrop, Scrop, Nt,
    %     reg_method);
    %     [RF{i},~,tmp] = SimpleGLM(Ecrop(i,:)', tcrop, Scrop, Nt, reg_method);
        [RF{i},~,tmp] = SimpleGLM(Ecrop(i,:)', tcrop, Scrop, regressor_range, reg_method);
        NewData{i} = tmp';
        % NewData(i,isRelevantTime) = tmp';
        fprintf('Done! took %.2f seconds.\n',toc)
    end

    for i=1:EEG.nbchan
        responseFns(:,:,i) = RF{i};
    end
end
% chansToDo = [16:20, 22, 25:29, 35, 56:61, 64, 67:76];
% chansToDo = [1:15, 21, 23:24, 30:34, 36:55, 62:63, 65:66, 77:79];
% parfor i=1:numel(chansToDo)
%     tic;
%     fprintf('Electrode %s...',chanLabels{chansToDo(i)})
%     [responseFns(:,:,i),~,tmp] = SimpleGLM(Ecrop(chansToDo(i),:)', tcrop, Scrop, Nt);
%     NewData{i} = tmp';
%     %     NewData(i,isRelevantTime) = tmp';
%     fprintf('Done! took %.2f seconds.\n',toc)
% end
% responseFns(:,:,chansToDo) = responseFns(:,:,1:numel(chansToDo));
% NewData(chansToDo) = NewData(1:numel(chansToDo));

responseFns = permute(responseFns,[3 2 1]); % electrodes x time x regressors
tResponse = (regressor_range(1):regressor_range(2))*dt;
NewEEG = EEG;
try
    NewReshapedData = zeros(size(reshapeddata));
    NewReshapedData(:,isRelevantTime) = cat(1,NewData{:});
    NewEEG.data = reshape(NewReshapedData,size(EEG.data));
catch
    warning('DJ:RunEegGlm:matrixWrite','Failed to write data as matrix... writing as cell array instead!');
    NewEEG.data = NewData;
end

% Plot results
if plotResults
    % Get channel number to plot
    iChan = find(strcmp(chanToPlot,{EEG.chanlocs(:).labels})); 
    
    % plot responses of each regressor
    for i=1:Nr
        MakeTopoMovie(responseFns(:,:,i),tResponse,EEG.chanlocs,20);
        MakeFigureTitle(sprintf('EEG response to regressor %s',event_types{i}));
    end

    % Plot one electrode's data
    figure;      
    handles = nan(1,Nr+1);
    handles(1) = subplot(Nr+1,1,1);
    plot(t,[reshapeddata(iChan,:); NewReshapedData(iChan,:)]');
    legend('true EEG','Reconstructed EEG');
    ylabel(sprintf('EEG (electrode %s)',chanToPlot));
    for i=1:Nr
        handles(i+1) = subplot(Nr+1,1,i+1);
        plot(t,s(i,:));
        ylabel(sprintf('Regressor %s',event_types{i}))    
    end
    xlabel('time (same units as input t)')
    linkaxes(handles,'x');
end