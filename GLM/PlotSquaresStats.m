function PlotSquaresStats(subject,setup_rule,during_suffix,after_suffix,eventsToPlot)

% Plot squares GLM response functions and single-subject stats.
%
% PlotSquaresStats(subject,setup_rule,during_suffix,after_suffix,eventsToPlot)
% PlotSquaresStats(EEG_before,during,after,after_suffix,eventsToPlot)
%
% INPUTS:
% -subject is a scalar indicating the subject number of your saved data.
% -setup_rule is the string used as input to SetUpGlm when you ran your
% GLM.
% -during_suffix is the end of the results file name saved after the
% nuisance GLM.
% -after_suffix is the end of the results file name saved after the final
% GLM.
% -eventsToPlot is a vector indicating the indices of the events you want
% to plot with PlotGlmResponses.
% -EEG_before is an eeg file with the proper events 
%
% Created 1/20/12 by DJ for one-time use.
% Updated 5/10/12 by DJ - made into a function
% Updated 6/27/12 by DJ - cascade figures across screen
% Updated 7/26/12 by DJ - no multiple comparisons correction

%% SACCADE-9REG (subtract out saccade regressor first)
% Handle defaults
if nargin<2 
    setup_rule = 'Saccade-Type-v1pt5';
end
if nargin<3
    during_suffix ='Saccade';
end
if nargin<4
    after_suffix = setup_rule;
end
if nargin<5
    eventsToPlot = 1:5;
end

% Load Data
if isstruct(subject) % Parse data
    EEG_before = subject;
    during = setup_rule;
    after = during_suffix;
    subject = str2double(EEG_before.subject);
else % Load data 
    EEG_before = pop_loadset(sprintf('sq-%d-all-filtered-50Hz-noduds.set',subject));
    EEG_before = SetUpGlm(subject,EEG_before,setup_rule);
    if ~isempty(during_suffix)
        during = load(sprintf('sq-%d-GLMresults-%s',subject,during_suffix));
    else
        during = [];
    end
    after = load(sprintf('sq-%d-GLMresults-%s',subject,after_suffix));
end



if ~isempty(during)
           
    % Remove nuisance events
    if isempty(during.nuisance_events) 
        during.nuisance_events = during.regressor_events; 
    end
    EEG_during = SubtractOutGlmResponses(EEG_before,during.responseFns,during.during.influence,nuisance_events); % will this work with old data?
%     after = after0;
%     after.regressor_events = after0.regressor_events([1 3 2 4:end]);
%     after.responseFns = after0.responseFns(:,:,[1 3 2 4:end]);    

else  
    EEG_during = EEG_before;       
end

% Get & plot single-subject z scores
% [z,p] = GetGlmZscore(EEG_during,after,'fdr');
% MakeFigureTitle(sprintf('sq-%d-GLMresults-%s, FDR corrected',subject,after_suffix));
[z,p] = GetGlmZscore(EEG_during,after,'none');
MakeFigureTitle(sprintf('sq-%d-GLMresults-%s, no multcompare correction',subject,after_suffix));

if ~isempty(eventsToPlot)
    % Plot results
    PlotGlmResponses(after.regressor_events(eventsToPlot), after.responseFns(:,:,eventsToPlot), ...
        after.tResponse, after.NewEEG.chanlocs, [], [-6 6]);
    PlotGlmResponses(after.regressor_events(eventsToPlot), z(:,:,eventsToPlot), ...
        after.tResponse, after.NewEEG.chanlocs, [], [-8 8]);
    % Arrange topomovie figures
    fh=findobj(0,'type','figure');
    CascadeFigures(fh(length(eventsToPlot)*2:-1:1),length(eventsToPlot));
end
