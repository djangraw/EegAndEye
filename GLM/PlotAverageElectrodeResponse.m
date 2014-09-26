function meaneegdata = PlotAverageElectrodeResponse(subjects,glmtype,electrode,eventnums,iLevel)

% Plot the response function for one electrode across event types.
%
% meaneegdata = PlotAverageElectrodeResponse(subjects,glmtype,electrode,eventnums,iLevel)
%
% INPUTS:
% -subjects is an n-element vector of subject numbers.
% -glmtype is a string indicating what kind of glm was done.  The data will
% be loaded from 'sq-<subjects(i)>-GLMresults-<glmtype>.mat', which must be
% in the current path.
% -electrode is a string indicating the electrode you want to see plotted.
% It must be present for all the subjects.
% -eventnums is a p-element vector of event numbers whose response
% functions will be plotted. (default: 1:5)
% -iLevel is a scalar indicating the level of analysis holding the event of
% interest (default = results.iLevel).
%
% OUTPUTS:
% -meaneegdata is a Txp matrix, where T is the number of time points in the
% response function, containing the mean response functions across subjects
% for each of the p events.
%
% Created 6/12/12 by DJ.
% Updated 6/18/12 by DJ - comments
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 3/19/14 by DJ - tResponse in cells

% Handle inputs
if nargin<4 || isempty(eventnums)
    eventnums = 1:5;
end
if nargin<5 
    iLevel = [];
end
nSubjects = numel(subjects);
nEvents = numel(eventnums);

% Load data
fprintf('Loading data')
for i=1:nSubjects
    fprintf('.')
    results = load(sprintf('sq-%d-GLMresults-%s',subjects(i),glmtype));
    a(i) = UpdateGlmResultsFormat(results);
    if isempty(iLevel)
        iLevel = a(i).iLevel;
    end
    % Convert to electrode number
    elecnum = find(strcmp(electrode,{a(i).NewEEG.chanlocs.labels}));
    
    eegdata(:,:,i) = squeeze(a(i).responseFns{iLevel}(elecnum,:,eventnums));
       
end
meaneegdata = mean(eegdata,3);
fprintf('done!\n')

disp('Plotting...')
% Plot data
cla; hold on;
% plot(a(1).tResponse,meaneegdata)
colors = {'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r--' 'g--' 'b--'};
if iscell(a(1).tResponse)
    tResponse = a(1).tResponse{a(1).iLevel};
else
    tResponse = a(1).tResponse;
end
for i=1:nEvents
    plot(tResponse,meaneegdata(:,i),colors{i})
end
    
% Annotate plot
legend(a(1).regressor_events{iLevel}(eventnums))
title(sprintf('subjects [%s], %s, %s',num2str(subjects),glmtype,electrode))
xlabel('time (ms)')
ylabel('response (uV)')
grid on
hold on
plot([tResponse(1), tResponse(end)],[0 0],'k')
plot([0 0],get(gca,'ylim'),'k')
disp('Done!')