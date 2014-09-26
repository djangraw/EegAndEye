function PlotScalpMaps(results,chanlocs,clim,tResults,rownames,cthresh)

% Plots a grid of scalp maps, where each row is a condition or subject and
% each column is a time window.
%
% PlotScalpMaps(results,chanlocs,clim,tResults,rownames,cthresh)
% 
% INPUTS:
% - results is a DxTxN matrix, where D = # channels T = # time windows, N =
%   # conditions.
% - chanlocs is a D-element vector of structs from an eeglab struct, i.e.,
%   EEG.chanlocs.
% - clim is a 2-element vector indicating the min and max color limits you
%   want to use on your plots. [default: [-max(abs(results(:))),
%   max(abs(results(:)))].
% - tResults is a T-element vector indicating the time of each time window
%   (it will be plotted above the plot as ("t=<time>"). [default: only put
%   row names in title]
% - rownames is an N-element cell array of strings indicating the title of
%   each row in your plot. [default: only put window times in title]
% - cthresh is a threshold below which colors will be set to the same color
%   as zero.
%
% Created 4/24/13 by DJ.
% Updated 4/29/13 by DJ - split into GetScalpMaps and this, added rownames 
%  input.
% Updated 5/23/13 by DJ - comments
% Updated 9/5/13 by DJ - one colorbar in lower left corner
% Updated 9/12/13 by DJ - added cthresh input
% Updated 8/13/14 by DJ - fixed clim=[0 0] bug
% Updated 8/18/14 by DJ - fixed clim=[-inf, inf] bug

% Handle defaults
if nargin<3 || isempty(clim)
    clim = [-max(abs(results(:))) max(abs(results(:)))];
end      
if nargin<4
    tResults = [];
end
if nargin<5 || isempty(rownames)    
    rownames = repmat({''},1,size(results,3));
end
if nargin<6 || isempty(cthresh)
    cthresh = 0;
end
if isequal(clim,[0 0])
    clim = [-1 1];
end
if isinf(clim(1)), clim(1) = -8; end
if isinf(clim(2)), clim(2) = 8; end
    

% Set up
[D,M,N] = size(results); % # elecs, # time points, # subjects/conditions
% make colormap
cmap = ThresholdedCmap(clim,cthresh);

TITLES_EVERY_PLOT = false;

for i=1:N % window
    for j=1:M % subject/condition
        % make topoplot
        subplot(N,M,(i-1)*M+j);
        topoplot(results(:,j,i),chanlocs,'plotrad',0.5,'electrodes','off');        
        % set plot colors
        set(gca,'clim',clim);
        colormap(cmap);
%         colorbar;
        if TITLES_EVERY_PLOT
            % add title
            if ~isempty(rownames{i}) && ~isempty(tResults) % plot both rownames and times
                title(sprintf('%s, t=%.1f ms',rownames{i},tResults(j)));
            elseif ~isempty(tResults) % just plot the times
                title(sprintf('t=%.1f ms',tResults(j)));
            elseif ~isempty(rownames{i}) % just plot the names
                title(rownames{i});            
            end % otherwise, no title!
        end
    end
end

if ~TITLES_EVERY_PLOT
    % Annotate
    if ~isempty(tResults)
        for i=1:M
            subplot(N,M,(N-1)*M+i);
            xlabel(sprintf('t = %.1f ms',tResults(i)),'visible','on');
        end
    end
    for j=1:N
        subplot(N,M,(j-1)*M+1);
        if ~isempty(rownames{j})
            ylabel(rownames{j},'visible','on');
        end
    end
end

% Put colorbar in lower right corner
lowerright = get(subplot(N,M,N*M),'Position');
axispos = [lowerright(1)+lowerright(3)/2, lowerright(2:4)];
axes('Position',axispos,'CLim',clim,'visible','off');%,'fontname',fontname,'fontsize',fontsize);
colormap(cmap);
colorbar%('fontname',fontname,'fontsize',fontsize)

    