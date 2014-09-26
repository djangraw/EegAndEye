function h = MakeTopoMovie(data,time,chanlocs,tbin_size,color_limits)

% Makes a figure/UI for scrolling through topoplot data like a movie.
%
% h = MakeTopoMovie(data,time,chanlocs,tbin_size)
%
% INPUTS:
%  data is a matrix of size [nElectrodes x nTimePoints] corresponding
% to the data we wish to plot at the times we wish to plot.
% - time is a nTimePoints-element vector containing the time in seconds
% at each point in data (e.g., the output of pop_comperp).
% - chanlocs is a struct taken from an EEGLAB struct's chanlocs field.
% - tbin_size is the size of time bins (in seconds). Default is no
% binning in time.
% - color_limits is a 2-element vector indicating the min and max values
% you want to define the colorbar scale.  Default is 
% [max(data(:)) -max(data(:))].
%
% OUTPUTS:
% - h is a struct containing the handles for various items on the
% figure.  It can be used to get or change properties of the figures.  For 
% example, type 'get(h.Plot)' to get the properties of the main topoplot.
%
% Created 8/9/10 by DJ.
% Updated 8/10/10 by DJ - reformatted to not be specific to subtracted ERPs
% Updated 8/12/10 by DJ - added input tbin_size, 'no binning' option, fixed
%                         time display formatting
% Updated 8/25/10 by DJ - now only 1 click on bottom plot will select time
% Updated 10/28/11 by DJ - added Save AVI button
% Updated 11/4/11 by DJ - fixed Save AVI frame rate bug
% Updated 12/29/11 by DJ - added color_limits input
% Updated 1/19/12 by DJ - fixed colorlim for data containing inf
% Updated 2/24/12 by DJ - comments

% Handle inputs
if nargin<4 || isempty(tbin_size), bin_data = 0; % no averaging over time
else bin_data = 1;
end
if nargin<5 || isempty(color_limits)
    color_max = max(max(data(~isinf(data))));
    color_limits = [-color_max color_max]; %topoplot color limits (in uV for an ERP, ERP ylabel displays these units)
end
% -------- DATA GATHERING -------- %

% Bin data
if bin_data
    disp('Binning data...')
    tbins = (time(1)+tbin_size/2) : tbin_size : (time(end)-tbin_size/2); % bins do not overlap
    erpbin = zeros(size(data,1),numel(tbins)); % preallocate matrix of erp.
    for i=1:numel(tbins)
        erpbin(:,i) = mean(data(:,time > tbins(i)-tbin_size/2 & time <= tbins(i)+tbin_size/2),2); % find data in the bin and average it
    end    
else
    tbins = time;
    erpbin = data;    
end
erpbin_original = erpbin; % make a static copy
color_thresholds = [0 0]; % colors within this range won't be plotted

% -------- INITIAL PLOTTING -------- %
disp('Setting up figure...');
figure; % make a new figure
iTime = 1; % the global index of the current time point - this is used throughout all functions

% Main scalp plot
h.Plot = axes('Units','normalized','Position',[0.15 0.25 0.8 0.7]); % set position
topoplot(erpbin(:,iTime),chanlocs,'MapLimits',color_limits,'conv','on'); % initial topoplot
colorbar; % add colorbar (this only needs to be called once)
title(sprintf('t = %.3f s',tbins(iTime)));

% ERP plot for selecting time
h.ERP = axes('Units','normalized','Position',[0.15 0.1 0.7 0.1]); % set position
plot(tbins,erpbin,'ButtonDownFcn',@erp_callback); % set up time selection function
xlabel('time (s)');
ylabel('ERP (uV)');
hold on
h.Line = plot(tbins([iTime iTime]), get(gca,'YLim'),'k','linewidth',2); % Line indicating current time
set(h.ERP,'ButtonDownFcn',@erp_callback) % this must be called after plotting, or it will be overwritten

% -------- GUI CONTROL SETUP -------- %
disp('Making GUI controls...')
h.Play = uicontrol('Style','togglebutton',...
                'String','Play',...
                'Units','normalized','Position',[.45 .2 .1 .05],...
                'Callback',@play_callback); % play button
h.Back = uicontrol('Style','pushbutton',...
                'String','<',...
                'Units','normalized','Position',[.38 .2 .05 .05],...
                'Callback',@back_callback); % back putton
h.Fwd = uicontrol('Style','pushbutton',...
                'String','>',...
                'Units','normalized','Position',[.57 .2 .05 .05],...
                'Callback',@fwd_callback); % forward button            
h.MPG = uicontrol('Style','pushbutton',...
                'String','Save AVI',...
                'Units','normalized','Position',[.15 .25 .15 .05],...
                'Callback',@saveavi_callback); % Make MPG button   
h.ColorLims = uicontrol('Style','slider',...
                'Min',0,...
                'Max',color_limits(end),...
                'Value',color_limits(end),...
                'Units','normalized','Position',[.05 .4 .03 .5],...
                'Callback',@color_callback); % color_limits slider
h.ColorThresh = uicontrol('Style','slider',...
                'Min',0,...
                'Max',color_limits(end),...
                'Value',0,...
                'Units','normalized','Position',[.15 .4 .03 .5],...
                'Callback',@color_callback); % color_thresholds slider
h.ColorLimLabel = uicontrol('Style','text',...
                'BackgroundColor',get(gcf,'Color'),...
                'String','CMax',...
                'Units','normalized','Position',[.02 .37 .09 .03]); % color_limits display
h.ColorThreshLabel = uicontrol('Style','text',...
                'BackgroundColor',get(gcf,'Color'),...
                'String','CThr',...
                'Units','normalized','Position',[.12 .37 .09 .03]); % color_limits display
h.ColorLimString = uicontrol('Style','edit',...
                'BackgroundColor',get(gcf,'Color'),...
                'String',sprintf('%.2f',color_limits(2)),...
                'Units','normalized','Position',[.02 .31 .09 .06],...
                'Callback',@color_callback); % color_limits display
h.ColorThreshString = uicontrol('Style','edit',...
                'BackgroundColor',get(gcf,'Color'),...
                'String',sprintf('%.2f',color_thresholds(2)),...
                'Units','normalized','Position',[.12 .31 .09 .06],...
                'Callback',@color_callback); % color_limits display
            
disp('Done!')
ginput(0); % For some reason, this solves the flashy problem!
    
% -------- SUBFUNCTIONS -------- %
function redraw() % Update the line and topoplot
    % Check that iTime is within allowable bounds
    if iTime<1 iTime=1;
    elseif iTime>size(erpbin,2) iTime = size(erpbin,2);
    end
    % Adjust plots
    set(h.Line,'XData',tbins([iTime iTime]));
    axes(h.Plot);
    topoplot(erpbin(:,iTime),chanlocs,'MapLimits',color_limits,'conv','on','electrodes','on');
    title(sprintf('t = %.3f s',tbins(iTime)));
end

function color_callback(hObject,eventdata)    
    if hObject==h.ColorLimString
        cmax = str2double(get(hObject,'String'));
        cmin = get(h.ColorThresh,'Value');
    elseif hObject==h.ColorThreshString
        cmin = str2double(get(hObject,'String'));
        cmax = get(h.ColorLims,'Value');
    else
        cmax = get(h.ColorLims,'Value');
        cmin = get(h.ColorThresh,'Value');
    end
    set(h.ColorLimString,'String',sprintf('%.2f',cmax));
    set(h.ColorThreshString,'String',sprintf('%.2f',cmin));
    set(h.ColorLims,'Value',cmax);
    set(h.ColorThresh,'Value',cmin);
    color_limits = [-cmax cmax];        
    color_thresholds = [-cmin cmin];
    erpbin = erpbin_original;
    erpbin(erpbin>color_thresholds(1) & erpbin<color_thresholds(2)) = 0;
    redraw();
end
    

function erp_callback(hObject,eventdata) % First mouse click on the ERP brings us here
    cp = get(h.ERP,'CurrentPoint'); % get the point(s) (x,y) where the person just clicked
    x = cp(1,1); % choose the x value of one point (the x values should all be the same).
%     [x,y] = ginput(1);
    iTime = find(tbins>=x,1); % find closest time to the click
    redraw; % update line and topoplot
end    

function play_callback(hObject,eventdata)
    % Get button value
    button_is_on = get(hObject,'Value') == get(hObject,'Max');
    % Keep incrementing and plotting
    while button_is_on && iTime < numel(tbins) %until we press pause or reach the end
        iTime=iTime+1;
        redraw;
%         drawnow;
        button_is_on = get(hObject,'Value') == get(hObject,'Max');
    end
    set(hObject,'Value',get(hObject,'Min')); % if we've reached the end, turn off the play button
end %function play_callback

function back_callback(hObject,eventdata)
    iTime = iTime-1; % decrement time index
    redraw; % update line and topoplot
end

function fwd_callback(hObject,eventdata)
    iTime = iTime+1; % increment time index
    redraw; % update line and topoplot
end

function saveavi_callback(hObject,eventdata)
    filename = uiputfile; % query user for filename
    if isequal(filename,0)
        return;
    end
    iTime = 0;
    fprintf(['Saving to ' filename '...'])
    aviobj = avifile(filename,'fps',0.1/(tbins(2)-tbins(1))); % slowed down 10x    
    while iTime<size(erpbin,2)
        if rem(iTime,100)==0
            fprintf('\n');
        end
        fprintf('.');
        fwd_callback;
        F = getframe(gcf);
        aviobj = addframe(aviobj,F);
    end
    fprintf('\nClosing...\n');
    aviobj = close(aviobj);
    disp('Done!')
end

end %function MakeTopoMovie