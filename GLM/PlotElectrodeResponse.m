function eegdata = PlotElectrodeResponse(subject,glmtype,electrode,eventnums,iLevel)

% Plot the response function for one electrode across event types.
%
% eegdata = PlotElectrodeResponse(subjects,glmtype,electrode,eventnums,iLevel)
%
% INPUTS:
% -subject is the subject number (scalar) or a loaded results struct.
% -glmtype is a string indicating what kind of glm was done.  The data will
% be loaded from 'sq-<subject>-GLMresults-<glmtype>.mat', which must be
% in the current path.
% -electrode is a string indicating the electrode you want to see plotted.
% It must be present for all the subjects.
% -eventnums is a p-element vector of event numbers whose response
% functions will be plotted.
% -iLevel is a scalar indicating the level of analysis holding the event of
% interest (default = results.iLevel).
%
% OUTPUTS:
% -eegdata is a Txp matrix, where T is the number of time points in the
% response function, containing the response functions for this subject
% for each of the p events.
%
% Created 6/12/12 by DJ.
% Updated 6/18/12 by DJ - comments
% Updated 8/10/12 by DJ - Update results format, NewEEG -> EEG
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 3/19/14 by DJ - tResponse in cells

% Handle inputs
if nargin<4 || isempty(eventnums)
    eventnums = 1:5;
end
if nargin<5
    iLevel = [];
end
nEvents = numel(eventnums);

% Load data
if isnumeric(subject)
    a = load(sprintf('sq-%d-GLMresults-%s',subject,glmtype));
    a = UpdateGlmResultsFormat(a);
else
    a = subject;
    a = UpdateGlmResultsFormat(a);
    subject = str2num(a.EEG.subject);
end
if isempty(iLevel)
    iLevel = a.iLevel;
end

% Convert to electrode number
if ischar(electrode)
    elecnum = find(strcmp(electrode,{a.EEG.chanlocs.labels}));
else
    elecnum = electrode;
end

% Plot data
eegdata = squeeze(a.responseFns{iLevel}(elecnum,:,eventnums));
cla; hold on;
% plot(a.tResponse,eegdata);
colors = {'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r--' 'g--' 'b--'};
if iscell(a.tResponse)
    tResponse = a.tResponse{a.iLevel};
else
    tResponse = a.tResponse;
end
for i=1:nEvents
    plot(tResponse,eegdata(:,i),colors{i})
end
ylim([-4 4])

% Annotate plot
legend(a.regressor_events{iLevel}(eventnums))
title(sprintf('S%d, %s, %s',subject,glmtype,a.EEG.chanlocs(elecnum).labels))
xlabel('time (ms)')
ylabel('response (uV)')
grid on
hold on
plot([tResponse(1), tResponse(end)],[0 0],'k')
plot([0 0],get(gca,'ylim'),'k')