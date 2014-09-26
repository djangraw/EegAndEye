function PlotSvdWeightsAndCourses(h,tH,weights,eigenvalues,chanlocs,reg_names,nComps,winv,colors)

% PlotSvdWeightsAndCourses(h,tH,weights,eigenvalues,chanlocs,reg_names,nComps,winv,colors))
%
% To be called by ApplySvdToGlmResults, or independently.
%
% INPUTS:
% -colors is an M-element cell array of strings or an Mx3 matrix indicating
%  the color for each line.
%
% Created 1/5/12 by DJ.
% Updated 7/18/12 by DJ - added nComps input
% Updated 1/31/13 by DJ - added winv input. topoplot winv, not weights
% Updated 3/13/14 by DJ - added Cmap input.

N = size(h,3);

if nargin<6 || isempty(reg_names)
    reg_names = cell(1,N);
    for i=1:N
        reg_names{i} = sprintf('Reg %d',i);
    end
end
if nargin<7 || isempty(nComps)
    nComps = 6;
end
if nargin<8 || isempty(winv)
    warning('No weight inverse given - plotting transpose of weights instead.')
    winv = weights'; 
end
if nargin<9 || isempty(colors)
%     colors = {'b' 'r' 'g' 'c' 'm' 'y' 'k' 'b--' 'r--' 'g--'};
    colors = distinguishable_colors(N,{'w','k'});
end
if ~iscell(colors)
    colorsmat = colors;
    colors = cell(1,size(colorsmat,1));
    for i=1:size(colorsmat,1)
        colors{i} = colorsmat(i,:);
    end
end


for i=1:min(nComps,length(chanlocs))
    % Plot topoplot
    subplot(nComps,2,2*i-1)
%     topoplot(weights(:,i),chanlocs,'electrodes','on');
    topoplot(winv(:,i),chanlocs,'electrodes','on');
    colorbar
    title(sprintf('component %d, lambda = %.4g',i,eigenvalues(i)));
    % Plot timecourse
    subplot(nComps,2,2*i)
    cla; hold on;
    for j=1:N
        plot(tH,weights(i,:)*h(:,:,j),'color',colors{j},'LineWidth',2);
    end
    set(gca,'xgrid','on')
    plot([tH(1) tH(end)],[0 0],'k-');
    plot([0 0],get(gca,'ylim'),'k-');
    ylabel('activity (uV)')
end
xlabel('time (ms)')
MakeLegend(colors(1:N),reg_names(1:N));

