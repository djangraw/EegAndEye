function PlotResponseFns(responseFns,regressor_events,tResponse,chanlocs,chanToPlot,colors)

% PlotResponseFns(responseFns,regressor_events,tResponse,chanlocs,chanToPlot,colors)
%
% INPUTS:
% -responseFn is a DxTxMxN matrix, where D=#chan, T=#samples, M=#events,
% and N=#subjects, containing the response function across all channels.
% -regressor_events is an M-element cell array of strings describing the
% events.
% -tResponse is a T-element vector containing the time of each sample
% relative to the locking event (in ms).
% -chanlocs is a D-element vector of structs from an eeglab file
% (EEG.chanlocs).
% -chanToPlot is a scalar indicating the channel number you want to plot,
% or a string indicating that channel's label in the chanlocs struct.
% -colors is an M-element cell array of strings or an Mx3 matrix indicating
% the color for each line.
%
% Created 2/7/13 by DJ.
% Updated 3/8/13 by DJ - added empty plot option to just get legend.
% Updated 7/16/13 by DJ - changed order to brgcmyk)
% Updated 10/8/13 by DJ - added colors input

if nargin<6 || isempty(colors)
    colors = {'b' 'r' 'g' 'c' 'm' 'y' 'k' 'b--' 'r--' 'g--'};
end
if ~iscell(colors)
    colorsmat = colors;
    colors = cell(1,size(colorsmat,1));
    for i=1:size(colorsmat,1)
        colors{i} = colorsmat(i,:);
    end
end


% Set up
nPlots = size(responseFns,3);
nSubjects = size(responseFns,4);
% Find index and name of electrode to plot
if ischar(chanToPlot)
    elecnum = find(strcmpi(chanToPlot,{chanlocs.labels}));
else
    elecnum = chanToPlot;
    chanToPlot = chanlocs(elecnum).labels;
end

% Plot data
eegdata = squeeze(mean(responseFns(elecnum,:,:,:),4));
eegstderr = squeeze(std(responseFns(elecnum,:,:,:),[],4))/sqrt(nSubjects);
if size(eegdata,1)==1, % row vectors mess up ErrorPatch
    eegdata = eegdata'; 
end


cla; hold on;
if ~isempty(eegdata)
    for i=1:nPlots
        plot(tResponse,eegdata(:,i),'color',colors{i},'LineWidth',2)
    end
    if nSubjects>1 % only do patches if there is a stderr
        for i=1:nPlots
            ErrorPatch(tResponse,eegdata(:,i)',eegstderr(:,i)',colors{i},colors{i});
        end
    end
else % Add dummy lines to get only the legend.
    for i=1:nPlots
        plot(tResponse(1)-100,0,'color',colors{i},'LineWidth',2)
    end
end
% ylim([-4 4])

% Annotate plot
legend(regressor_events)
title(sprintf('%d subjects, %s',nSubjects,chanToPlot))
xlabel('time (ms)')
ylabel('response +/- stderr (uV)')
grid on
hold on
plot([tResponse(1), tResponse(end)],[0 0],'k')
plot([0 0],get(gca,'ylim'),'k')
xlim([tResponse(1), tResponse(end)]);