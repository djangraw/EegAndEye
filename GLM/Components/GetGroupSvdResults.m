function [weights, timecourse, avgW, avgTC, chanlocs, tResponse, steW, steTC, eventnames] = GetGroupSvdResults(S,nComponents,iEventTypes,timeWindow,componentSwap,doPlot)

% Gets and plots the average +/- std dev of the weights and timecourses
% across subjects of the top 2 SVD components.
%
% [weights, timecourse, avgW, avgTC, chanlocs, tResponse, steW, steTC,
% eventnames] = GetGroupSvdResults(S,nComponents,iEventTypes,timeWindow,
% componentSwap)
% 
% INPUTS:
% - S is an n-element vector of GLM results from the same analysis on n
% different subjects.
% - nComponents is the number of top components you want to use.  The 1st
% component of all subjects will be considered equivalent.  Same with the
% others.
% - iEventTypes is a vector of the event type indices you want to use. If 
% it is a scalar, then events 1:iEventTypes will be used. The event type
% order should be the same across all subjects. 
% - timeWindow should be a 2-element vector indicating the start and end
% time (in ms) you want to use when calculating the SVD components.  See
% ApplySvdToGlmResults.m for details.
% - componentSwap is a matrix of size nComponents x n and specify the order
% and sign of the weights matrices to be used in the group results.  For
% example, componentSwap=[1 1; 2 -3] means that the 1st group component
% will be the average of both subjects' 1st components, but the 2nd group
% component will be the average of S1's 2nd component and the inverse of 
% S2's 3rd component.
%
% OUTPUTS:
% - weights is an n-element cell array of nChannels x nComponents matrices.
%  Each column is a set of weights for that subject and component.
% - timecourse is an nComponents x nEventTypes cell array of nSubjects x T
% matrices (where T is the number of time points in timeWindow).  Each row
% is the time course of activity for that subject, component, and event
% type.
% - avgW is the average weights over all the subjects
% - avgTC is the average timecourse over all the subjects
% - chanlocs is a 1xnChan struct array of channel locations that can be used as
% input to topoplot, etc.
% - tResponse is a T-element vector of the times corresponding to each
% timecourse element.
% - steW and steTC are the standard errors of the weights and timecourses.
%
% Created 4/17/12 by DJ.
% Updated 4/23/12 by DJ - updated component sign changes
% Updated 5/8/12 by DJ - added outputs, while loop in weight adjustment
% Updated 5/9/12 by DJ - changed fontsize for SfN abstract plotting
% Updated 5/11 12 by DJ - added componentSwap input
% Updated 7/18/12 by DJ - added nComps input
% Updated 7/20/12 by DJ - made compatible with new RunGlmGui results
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 10/2/13 by DJ - added tResponse, eventnames, stderr outputs for 
%   easy input to plot fn
% Updated 10/10/13 by DJ - changed nEventTypes input to iEventTypes
% Updated 3/19/14 by DJ - tResponse in cells

nSubjects = numel(S);
if nargin<2 || isempty(nComponents)
    nComponents = 3;
end
if nargin<3 || isempty(iEventTypes)
    iEventTypes = 1:5;    
end
if nargin<4 || isempty(timeWindow)
    timeWindow = []; % [200 500]
end
if nargin<5 || isempty(componentSwap)
    componentSwap = ones(nComponents,nSubjects);
    for i=1:nComponents
        componentSwap(i,:) = i;
    end
end
if nargin<6 || isempty(doPlot)
    doPlot = 1;
end
% Allow scalar input for iEventTypes
if numel(iEventTypes) == 1
    nEventTypes = iEventTypes;
    iEventTypes = 1:nEventTypes;
else
    nEventTypes = numel(iEventTypes);
end

% Update datasets
for i=1:numel(S)
    S_new(i) = UpdateGlmResultsFormat(S(i));
end
S = S_new;
clear S_new

%% Get weights
weights = ApplySvdToMultiSubjects([],S,timeWindow,iEventTypes,nComponents,doPlot);
% Make specified changes to component weights
weights0 = weights;
for i=1:nSubjects
    for j=1:nComponents
        weights{i}(:,j) = sign(componentSwap(j,i))*weights0{i}(:,abs(componentSwap(j,i)));
    end
%     weights{i}(:,2) = -weights{i}(:,2); % TEMP!
end

% Adjust so all components have the "same sign"
isok = zeros(nSubjects,nComponents);
while ~all(isok(:))
    isok(:) = 1;
    [avgW_temp,~,chanlocs_temp] = GetAverageWeights([S.EEG],weights);
    for i=1:nSubjects
        iChans = ismember({S(i).EEG.chanlocs.labels},{chanlocs_temp.labels}); % find channels included in group
        for j=1:nComponents
%             fprintf('S%s, c%d: %g\n',S(i).EEG.subject,j,weights{i}(iChans,j)'*avgW_temp(:,j));
            if weights{i}(iChans,j)'*avgW_temp(:,j)<0 % check whether individual and group are coherent
                fprintf('Reversing sign on subject %s, component %d...\n',S(i).EEG.subject,j);
                weights{i}(:,j) = -weights{i}(:,j); % if not, reverse sign on individual's weights
                isok(i,j) = 0;
            end
        end
    end
end

%% Get timecourses
timecourse = cell(nComponents,nEventTypes);
for i=1:nSubjects % subject
    for j=1:nComponents % component
        for k=1:nEventTypes % trial type
            timecourse{j,k}(i,:) = weights{i}(:,j)'*S(i).responseFns{S(i).iLevel}(:,:,iEventTypes(k));
        end
    end
end

%% Get average component
[avgW,steW,chanlocs] = GetAverageWeights([S.EEG],weights);

%% Get averages and standard errors
avgTC = zeros(nEventTypes,size(timecourse{1},2),nComponents);
steTC = zeros(size(avgTC));
for j=1:nComponents
    for k=1:nEventTypes
        avgTC(k,:,j) = mean(timecourse{j,k},1);
        steTC(k,:,j) = std(timecourse{j,k},[],1)/sqrt(nSubjects);
    end
end

%% Plot average and standard errors
tResponse = S(1).tResponse{S(1).iLevel};
eventnames = S(1).regressor_events{S(1).iLevel}(iEventTypes);
if doPlot
    figure;
    PlotGroupSvdResults(avgW, avgTC, chanlocs, tResponse, steW, steTC, eventnames);
end
% figure(877); clf;
% colors = {'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r--' 'g--' 'b--'};
% figure(877); clf;
% for j=1:nComponents
%     % plot weights
%     subplot(nComponents,4,j*4-3);
%     set(gca,'FontSize',15);
%     topoplot(avgW(:,j),chanlocs,'electrodes','on');
%     title(sprintf('Component %d: Mean Weights',j))
%     set(gca,'CLim',[-.2 .2])
%     colorbar('FontSize',15)
%     % plot stderr of weights
%     subplot(nComponents,4,j*4-2);
%     set(gca,'FontSize',15);
%     topoplot(steW(:,j),chanlocs,'electrodes','on');
%     title(sprintf('Component %d: StdErr Weights',j))
%     set(gca,'CLim',[-.1 .1])
%     colorbar('FontSize',15)
%     % plot timecourse
%     subplot(nComponents,2,j*2); cla; hold on;
%     for k=1:nEventTypes
%         plot(S(1).tResponse,avgTC(k,:,j),colors{k})
%     end
%     legend(S(1).regressor_events{S(1).iLevel}(1:nEventTypes),'Location','NorthWest')
%     for k=1:nEventTypes
%         ErrorPatch(S(1).tResponse,avgTC(k,:,j),steTC(k,:,j),colors{k},colors{k});
%     end
%     % annotate plot
%     set(gca,'xgrid','on','box','on','FontSize',15)
%     plot([S(1).tResponse(1) S(1).tResponse(end)],[0 0],'k-');
%     plot([0 0],get(gca,'ylim'),'k-');
%     ylabel('activity (uV)')
%     xlabel('time relative to saccade landing (ms)')
%     title(sprintf('Component %d: Mean & StdErr Timecourse',j))
% end
% % MakeLegend(colors([1 3 2 4]),{'Low Anticipation Non-Target', ...
% %     'Low Anticipation Target','High Anticipation Non-Target','High Anticipation Target'},[2 2 2 2])
% set(gca,'FontSize',15)

function [avgW, steW, chanlocs] = GetAverageWeights(ALLEEG,weights)
% Set up
nSubjects = length(weights);
nComponents = min([ALLEEG.nbchan]);
% Get list of common channels
channels = {ALLEEG(1).chanlocs.labels};
for i=1:nSubjects
    channels = intersect(channels, {ALLEEG(i).chanlocs.labels});
end
chanlocs = ALLEEG(1).chanlocs(ismember({ALLEEG(1).chanlocs.labels},channels));
% Get average weights
wts = cell(1,nComponents);
for i=1:nSubjects
    iChans = ismember({ALLEEG(i).chanlocs.labels},channels);
    for j=1:nComponents
        wts{j}(:,i) = weights{i}(iChans,j);
    end
end
avgW = zeros(numel(channels),nComponents);
steW = zeros(size(avgW));
for j=1:nComponents
    avgW(:,j) = mean(wts{j},2);
    steW(:,j) = std(wts{j},[],2)/sqrt(nSubjects);
end
