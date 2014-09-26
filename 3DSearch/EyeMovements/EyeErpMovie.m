function hOut = PlotEyeErps(data)

% Makes a figure/UI for scrolling through eye position data.
%
% PlotEyeErps(eyepos,pupilsize,x)
%
% The circle represents the subject's final eye position.  The dot's size represents
% the reported pupil size.  The big rectangle is the limits of the screen.
% Unity reports objects that are off the screen, but they aren't visible
% to the subject until they enter this rectangle.  The smaller rectangles
% represent the screen bounds of objects in the scene, as reported by a
% Unity replay.
%
% Inputs:
%   - eyepos is an nx2 matrix, where n is the number of samples of eye
% position data. Each row is the (x,y) position of the eye at that time.
% This will be the position of the dot on the screen.
%   - pupilsize is an n-element vector, where each element is the size of
% the subject's pupil (in unknown units).  This will be the size of the dot
% on the screen.
%   - x is a 3DSearch data structure as imported by Import_3DS_Data_v3.
%
% Outputs:
%   - hOut is a struct containing the handles for various items on the
% figure.  It can be used to get or change properties of the figures.  For 
% example, type 'get(h.Plot)' to get the properties of the main movie.
%
% Created 3/28/11 by DJ.
% TO DO: Comment, combine multiple sessions.

% %% -------- SETUP -------- %
% fs = 1000;
% screen_res = [1024, 768];
% % global epochSamples h;
% epochSamples = -500:1000;
% epochTimes = epochSamples*1000/fs;
% % Parse inputs
% object_limits = x.eyelink.object_limits;
% record_time = x.eyelink.record_time; % time 'offset' - the time when eyelink started recording.
% isTarget = strcmp('TargetObject',{x.objects(:).tag});
% 
% % -------- EPOCH -------- %
% % Get sample numbers around when targets or distractors came onscreen
% GetNumbers;
% targetOnset = x.eyelink.object_events(ismember(x.eyelink.object_events(:,2),Numbers.ENTERS+find(isTarget)),1);
% distractorOnset = x.eyelink.object_events(ismember(x.eyelink.object_events(:,2),Numbers.ENTERS+find(~isTarget)),1);
% for i=1:numel(targetOnset)
%     targetSamples(i,:) = targetOnset(i) + epochSamples;
% end
% for i=1:numel(distractorOnset)
%     distractorSamples(i,:) = distractorOnset(i) + epochSamples;
% end
% 
% targetEyeEpochs = cell(size(targetSamples));
% targetPupEpochs = zeros(size(targetSamples));
% targetObjEpochs = cell(size(targetSamples));
% distractorEyeEpochs = cell(size(distractorSamples));
% distractorPupEpochs = zeros(size(distractorSamples));
% distractorObjEpochs = cell(size(distractorSamples));
% 
% for i=1:numel(targetSamples)
%     targetEyeEpochs{i} = eyepos(targetSamples(i)-record_time+1,:)';
%     targetPupEpochs(i) = pupilsize(targetSamples(i)-record_time+1);
%     iLim = find(object_limits(:,1)>targetSamples(i)-20 & object_limits(:,1)<targetSamples(i),1,'last');
%     targetObjEpochs{i} = object_limits(iLim,3:6);
% end
% for i=1:numel(distractorSamples)
%     distractorEyeEpochs{i} = eyepos(distractorSamples(i)-record_time+1,:)';
%     distractorPupEpochs(i) = pupilsize(distractorSamples(i)-record_time+1);
%     iLim = find(object_limits(:,1)>distractorSamples(i)-20 & object_limits(:,1)<distractorSamples(i),1,'last');
%     distractorObjEpochs{i} = object_limits(iLim,3:6);
% end
% 
% % add to data struct
% data.epochSamples = epochSamples;
% data.epochTimes = epochTimes;
% data.targetEyeEpochs = targetEyeEpochs;
% data.targetPupEpochs = targetPupEpochs;
% data.targetObjEpochs = targetObjEpochs;
% data.distractorEyeEpochs = distractorEyeEpochs;
% data.distractorPupEpochs = distractorPupEpochs;
% data.distractorObjEpochs = distractorObjEpochs;
% data.screen_res = screen_res;

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
h.TargPlot = axes('Units','normalized','Position',[0.13 0.3 0.35 0.65]); % set position
rectangle('Position',[0 0 screen_res]);
title('Target Trials');
h.DistPlot = axes('Units','normalized','Position',[0.53 0.3 0.35 0.65]); % set position
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
    tEye = h.data.targetEyeEpochs(:,(iStart:iEnd)-epochSamples(1));
    tObj = h.data.targetObjEpochs(:,(iStart:iEnd)-epochSamples(1));
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
        if ~isempty(trialObjEnd)
            rectangle('Position',[trialObjEnd],'EdgeColor','k');
        end
        % plot eye position
        trialEye = [tEye{i,:}];
        plot(trialEye(1,:),trialEye(2,:),'b-');
        plot(trialEye(1,end),trialEye(2,end),'bo');
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
        plot(trialEye(1,:),trialEye(2,:),'r-');
        plot(trialEye(1,end),trialEye(2,end),'ro');
    end
    rectangle('Position',[0 0 screen_res]);
end


function new_position = roi_constraint(current_position) % First mouse click on the Time plot brings us here
    % constrain position
    h = get(gcf,'UserData');
    new_position = [current_position(1), min(get(h.Time,'YLim')), current_position(3), range(get(h.Time,'YLim'))];
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





