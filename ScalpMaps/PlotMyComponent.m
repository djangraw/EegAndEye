function [trialActivation, component] = PlotMyComponent(component,EEG,trials)

% Plot component map of given component or time ERP.
%
% trialActivation = PlotMyComponent(component,EEG,trials)
% [trialActivation, component] = PlotMyComponent(iSample,EEG,trials)
%
% INPUTS:
% - component is a Dx1 vector of electrode weights defining the component.
% - EEG is an eeglab data struct, in which EEG.data is a D x T x n matrix
% of voltages (mV).
% - trials is a vector of the indices of trials that should be included in
% the component ERP.
% - iSample is a scalar indicating the sample at which the ERP of all 
% trials should be considered the component.
%
% OUTPUTS:
% - trialActivation is an n x T matrix containing the activation of the
% given component at a given time.
% - component is a Dx1 vector output in case the iSample input option is
% used.
%
% Created 3/30/11 by DJ.
% Updated 3/31/11 by DJ - gives component as output, plots button presses,
% has color limits
% Updated 4/14/11 by DJ - added t=0 and V=0 lines
% Updated 9/23/14 by DJ - comments.

% set parameters
plotButton = 1; % plot button press times as black dots
colorlimits = [-200 200]; % set to empty for automatic

% Handle inputs
if nargin<3 || isempty(trials)
    nTrials = size(EEG.data,3);
    trials = 1:nTrials;
else
    nTrials = numel(trials);
end

if plotButton
    buttonTime = NaN(1,nTrials);
    GetNumbers;
end

if numel(component)==1   
    iSample = component;
    % Get component from ERP at certain sample number
    component = mean(EEG.data(:,iSample,:),3);
    titlestring = sprintf('component from ERP, t=%.1f ms',EEG.times(iSample));
else
    titlestring = 'given component';
end

% normalize component
component = component / sqrt(sum(component.^2)); % squared values should sum to 1

% plot component
subplot(2,2,1)
topoplot(component,EEG.chanlocs);
title(titlestring);

% get single-trial activations
trialActivation = zeros(nTrials,size(EEG.data,2));
for i=1:nTrials
    trialActivation(i,:) = component'*EEG.data(:,:,trials(i));
    if plotButton
        if sum(strcmp(EEG.epoch(trials(i)).eventtype,num2str(Numbers.BUTTON)))==1 % if there was a button press
            buttonTime(i) = EEG.epoch(trials(i)).eventlatency{strcmp(EEG.epoch(trials(i)).eventtype,num2str(Numbers.BUTTON))};
        end
    end
end

% plot activation ERP
subplot(2,2,2)
hold on
plot(EEG.times,mean(trialActivation,1));
xlabel('time (ms)');
ylabel('voltage (mV)');
title('Average component activation')
% plot lines at time t=0 and y=0
plot([0,0],get(gca,'YLim'),'k');
plot(get(gca,'XLim'),[0,0],'k:');

% plot single-trial activations
subplot(2,1,2)
imagesc(EEG.times,1:nTrials,trialActivation);
colorbar
xlabel('time (ms)');
ylabel('trial number');
title('Trial by trial component activation')
if ~isempty(colorlimits)
    set(gca,'CLim',colorlimits);
end
% plot line at time t=0
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');
hold on
plot([0,0],get(gca,'YLim'),'k');
set(gca,'XLim',xlims,'YLim',ylims);

% plot button press times
if plotButton
    hold on
    plot(buttonTime,1:nTrials,'k.')
end