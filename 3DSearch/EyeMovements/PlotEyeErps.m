function hOut = PlotEyeErps(data)

% Makes a figure/UI for scrolling through eye position data.
%
% hOut = PlotEyeErps(data)
%
% The circle represents the subject's final eye position.  The dot's size represents
% the reported pupil size.  The big rectangle is the limits of the screen.
% Unity reports objects that are off the screen, but they aren't visible
% to the subject until they enter this rectangle.  The smaller rectangles
% represent the screen bounds of objects in the scene, as reported by a
% Unity replay.
%
% Inputs:
%   - data is a struct from GetEyeErps.m containing info about times, eye
% position, etc. (see GetEyeErps for details)
%
% Outputs:
%   - hOut is a struct containing the handles for various items on the
% figure.  It can be used to get or change properties of the figures.  For 
% example, type 'get(h.Plot)' to get the properties of the main movie.
%
% Created 3/28/11 by DJ.
% Updated 4/1/11 by DJ - get input from GetEyeErps.m
% Updated 5/9/13 by DJ - flipped y axis (0 at the top!)


% -------- UNPACK -------- %
epochSamples = data.epochSamples;
screen_res = data.screen_res;

% -------- INITIAL PLOTTING -------- %
disp('Setting up figure...');
figure; % make a new figure
% global iStart iEnd; % the global index of the current time point - this is used throughout all functions
iStart = 1; 
iEnd = 2;

% Main eye plots
h.TargPlot = axes('Units','normalized','Position',[0.13 0.3 0.35 0.65],'ydir','reverse'); % set position
rectangle('Position',[0 0 screen_res]);
title('Target Trials');
h.DistPlot = axes('Units','normalized','Position',[0.53 0.3 0.35 0.65],'ydir','reverse'); % set position
rectangle('Position',[0 0 screen_res]);
title('Distractor Trials');
axis(h.TargPlot,[-200 screen_res(1)+200 -200 screen_res(2)+200]);
axis(h.DistPlot,[-200 screen_res(1)+200 -200 screen_res(2)+200]);

% title(sprintf('t = %.3f s',epochTimes(iTime)));

% -------- TIME SELECTION PLOT SETUP -------- %
% 'Time plot' for selecting and observing current time
h.Time = axes('Units','normalized','Position',[0.13 0.1 0.775 0.1],'Yticklabel','','ButtonDownFcn',@time_callback); % set position
% set(h.Time,'Xticklabel',get(gca,'XTick'));
hold on;
% Annotate plot
xlim([epochSamples(1) epochSamples(end)]);
xlabel('time (samples)');
% h.Rect = rectangle('Position',[iStart,
% min(get(h.Time,'YLim')),iEnd-iStart,range(get(h.Time,'YLim'))],'facecolor','c','ButtonDownFcn',@time_callback); % Line indicating current time
plot([0 0],get(gca,'YLim'),'k','ButtonDownFcn',@time_callback);
h.Rect = imrect(h.Time,[iStart, min(get(h.Time,'YLim')),iEnd-iStart,range(get(h.Time,'YLim'))],...
    'PositionConstraintFcn',@roi_constraint);
addNewPositionCallback(h.Rect,@roi_callback); % Line indicating current time
set(h.Time,'ButtonDownFcn',@time_callback); % this must be called after plotting, or it will be overwritten
title('Drag rectangle to set limits of plot')

% -------- GUI CONTROL SETUP -------- %
% disp('Making GUI controls...')
% h.Play = uicontrol('Style','togglebutton',...
%                 'String','Play',...
%                 'Units','normalized','Position',[.45 .2 .1 .05],...
%                 'Callback',@play_callback); % play button
% h.Speed = uicontrol('Style','slider',...
%                 'Min',1,'Max',100,'Value',1,'SliderStep',[.05 .2],...
%                 'Units','normalized','Position',[.45 .25 .1 .025]); % speed slider
% h.Back = uicontrol('Style','pushbutton',...
%                 'String','<',...
%                 'Units','normalized','Position',[.38 .2 .05 .05],...
%                 'Callback',@back_callback); % back putton
% h.Fwd = uicontrol('Style','pushbutton',...
%                 'String','>',...
%                 'Units','normalized','Position',[.57 .2 .05 .05],...
%                 'Callback',@fwd_callback); % forward button            
     
h.data = data;
set(gcf,'UserData',h);
hOut = h;
disp('Done!')
end %function MakeEyeMovie
    
% -------- SUBFUNCTIONS -------- %
function redraw() % Update the lines and rectangles

    % unpack data
    h = get(gcf,'UserData');
    iStart = h.data.iStart;
    iEnd = h.data.iEnd;
    epochSamples = h.data.epochSamples;
    screen_res = h.data.screen_res;
    
    % Check that iTime is within allowable bounds
    if iStart<epochSamples(1), iStart = epochSamples(1); end
    if iEnd>epochSamples(end), iEnd = epochSamples(end); end
    if iEnd<=iStart, iEnd = iStart+1; end
    
    % unpack and crop more data
    tSess = h.data.targetEventSessions;
    tTime = h.data.targetEventTimes;
    tEye = h.data.targetEyeEpochs(:,(iStart:iEnd)-epochSamples(1));
    tObj = h.data.targetObjEpochs(:,(iStart:iEnd)-epochSamples(1));
    dSess = h.data.distractorEventSessions;
    dTime = h.data.distractorEventTimes;
    dEye = h.data.distractorEyeEpochs(:,(iStart:iEnd)-epochSamples(1));
    dObj = h.data.distractorObjEpochs(:,(iStart:iEnd)-epochSamples(1));
    
    % Adjust plots
    title(h.Time,sprintf('From t=%.1f to t=%.1f ms',h.data.epochTimes(epochSamples==iStart),h.data.epochTimes(epochSamples==iEnd)));
%     set(h.Rect,'Position',[iStart, min(get(h.Time,'YLim')),iEnd-iStart,range(get(h.Time,'YLim'))]);
%     plot([0 0],get(gca,'YLim'),'k','ButtonDownFcn',@time_callback);
    
    % Plot lines for target and distractor trials
    axes(h.TargPlot); cla; hold on;   
    for i=1:size(tEye,1) % for each target trial
        % plot object position
        trialObjEnd = tObj{i,end};
        if ~isempty(trialObjEnd) && all(trialObjEnd(3:4)>0)
            rectangle('Position',[trialObjEnd],'EdgeColor','k');
        end
        % plot eye position
        trialEye = [tEye{i,:}];
%         plot(trialEye(1,:),trialEye(2,:),'b-');
%         plot(trialEye(1,end),trialEye(2,end),'bo');
        plot(trialEye(1,:),trialEye(2,:),'b-','UserData',...
            sprintf('clicked target trial %d (session %d, time %.3f)',i,tSess(i),tTime(i)),...
            'ButtonDownFcn','disp(get(gco,''UserData''))');
        plot(trialEye(1,end),trialEye(2,end),'bo','UserData',...
            sprintf('clicked target trial %d (session %d, time %.3f)',i,tSess(i),tTime(i)),...
            'ButtonDownFcn','disp(get(gco,''UserData''))');
    end
    rectangle('Position',[0 0 screen_res]);
    
    axes(h.DistPlot); cla; hold on;
    for i=1:size(dEye,1) % for each distractor trial
        % plot object position
        trialObjEnd = dObj{i,end};
        if ~isempty(trialObjEnd)
            rectangle('Position',[trialObjEnd],'EdgeColor','k');
        end
        % plot eye position
        trialEye = [dEye{i,:}];
%         plot(trialEye(1,:),trialEye(2,:),'r-');
%         plot(trialEye(1,end),trialEye(2,end),'ro');
        plot(trialEye(1,:),trialEye(2,:),'r-','UserData',...
            sprintf('clicked distractor trial %d (session %d, time %.3f)',i,dSess(i),dTime(i)),...
            'ButtonDownFcn','disp(get(gco,''UserData''))');
        plot(trialEye(1,end),trialEye(2,end),'ro','UserData',...
            sprintf('clicked distractor trial %d (session %d, time %.3f)',i,dSess(i),dTime(i)),...
            'ButtonDownFcn','disp(get(gco,''UserData''))');
    end
    rectangle('Position',[0 0 screen_res]);
end


function new_position = roi_constraint(current_position) % First mouse click on the Time plot brings us here
    % constrain position
    h = get(gcf,'UserData');
    xStart = max(current_position(1), min(get(h.Time,'XLim'))+1); % not too early
    xStart = min(xStart, max(get(h.Time,'XLim'))-current_position(3)); % not too late
    new_position = [xStart, min(get(h.Time,'YLim')), current_position(3), range(get(h.Time,'YLim'))];
end

function roi_callback(current_position)
    % get time points
    iStart = round(current_position(1));
    iEnd = round(current_position(1)+current_position(3));    
    % Update userdata
    h = get(gcf,'UserData');
    h.data.iStart = iStart;
    h.data.iEnd = iEnd;
    set(gcf,'UserData',h);
    % update lines and rectangles
    redraw; 
end   

function time_callback(hObject,eventdata) % First mouse click on the Time plot brings us here
%     % get new rectangle
%     disp('Click and drag rectangle across time points of interest')
%     rect = getrect(gca);
%     iStart = round(rect(1));
%     iEnd = round(rect(1)+rect(3));
%     % Update userdata
%     h = get(gcf,'UserData');
%     h.data.iStart = iStart;
%     h.data.iEnd = iEnd;
%     set(gcf,'UserData',h);
%     % update lines and rectangles
%     redraw; 
end    

% function play_callback(hObject,eventdata)
%     % Get button value
%     button_is_on = get(hObject,'Value') == get(hObject,'Max');
%     % Keep incrementing and plotting
%     while button_is_on && iTime < numel(time) %until we press pause or reach the end
%         iTime=iTime+floor(get(h.Speed,'Value'));
%         redraw;
%         button_is_on = get(hObject,'Value') == get(hObject,'Max');
%     end
%     set(hObject,'Value',get(hObject,'Min')); % if we've reached the end, turn off the play button
% end %function play_callback
% 
% function back_callback(hObject,eventdata)
%     iTime = iTime-1; % decrement time index
%     redraw; % update line and topoplot
% end
% 
% function fwd_callback(hObject,eventdata)
%     iTime = iTime+1; % increment time index
%     redraw; % update line and topoplot
% end





