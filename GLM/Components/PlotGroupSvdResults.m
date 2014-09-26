function PlotGroupSvdResults(avgW, avgTC, chanlocs, tResponse, steW, steTC, eventnames,Cmap)

% PlotGroupSvdResults(avgW, avgTC, chanlocs, tResponse, steW, steTC, eventnames,Cmap)
%
% Created 10/2/13 by DJ.

if nargin<4 || isempty(steW)
    steW = zeros(size(avgW));
end
if nargin<5 || isempty(steTC)
    steTC = zeros(size(avgTC));
end

nEventTypes = size(avgTC,1);
nComponents = size(avgTC,3);
if nargin<8 || isempty(Cmap)
    Cmap = distinguishable_colors(nEventTypes,{'w','k'}); % don't allow black or white
end
% colors = {'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r--' 'g--' 'b--'};
% figure(877); clf;
clf;
for j=1:nComponents
    % plot weights
    subplot(nComponents,4,j*4-3);
    set(gca,'FontSize',15);
    topoplot(avgW(:,j),chanlocs,'electrodes','on');
    title(sprintf('Component %d: Mean Weights',j))
    set(gca,'CLim',[-.2 .2])
    colorbar('FontSize',15)
    % plot stderr of weights
    subplot(nComponents,4,j*4-2);
    set(gca,'FontSize',15);
    topoplot(steW(:,j),chanlocs,'electrodes','on');
    title(sprintf('Component %d: StdErr Weights',j))
    set(gca,'CLim',[-.1 .1])
    colorbar('FontSize',15)
    % plot timecourse
    subplot(nComponents,2,j*2); cla; hold on;
    for k=1:nEventTypes
        plot(tResponse,avgTC(k,:,j),'Color',Cmap(k,:))
    end
%     legend(eventnames,'Location','NorthWest')
    for k=1:nEventTypes
        ErrorPatch(tResponse,avgTC(k,:,j),steTC(k,:,j),Cmap(k,:),Cmap(k,:));
    end
    % annotate plot
    set(gca,'xgrid','on','box','on','FontSize',15)
    plot([tResponse(1) tResponse(end)],[0 0],'k-');
    plot([0 0],get(gca,'ylim'),'k-');
    ylabel('activity (uV)')
    xlabel('time relative to saccade landing (ms)')
    title(sprintf('Component %d: Mean & StdErr Timecourse',j))
end
if isnumeric(Cmap)
    Cmap_cell = cell(1,size(Cmap,1));
    for i=1:size(Cmap,1)
        Cmap_cell{i} = Cmap(i,:);
    end
    MakeLegend(Cmap_cell,eventnames,2);
else
    MakeLegend(Cmap,eventnames,2);
end
% MakeLegend(colors([1 3 2 4]),{'Low Anticipation Non-Target', ...
%     'Low Anticipation Target','High Anticipation Non-Target','High Anticipation Target'},[2 2 2 2])
set(gca,'FontSize',15)