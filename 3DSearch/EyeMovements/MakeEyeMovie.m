function h = MakeEyeMovie(samples,pupilsize,x,t_start)

% Makes a figure/UI for scrolling through eye position data like a movie.
%
% MakeEyeMovie(samples,pupilsize,x)
% MakeEyeMovie(subject,session)
%
% The dot represents the subject's eye position.  The dot's size represents
% the reported pupil size.  The big rectangle is the limits of the screen.
% Unity reports objects that are off the screen, but they aren't visible
% to the subject until they enter this rectangle.  The smaller rectangles
% represent the screen bounds of objects in the scene, as reported by a
% Unity replay.
%
% Inputs:
%   - samples is an nx2 matrix, where n is the number of samples of eye
% position data. Each row is the (x,y) position of the eye at that time.
% This will be the position of the dot on the screen.
%   - pupilsize is an n-element vector, where each element is the size of
% the subject's pupil (in unknown units).  This will be the size of the dot
% on the screen.
%   - x is a 3DSearch data structure as imported by Import_3DS_Data_v2.
%   - subject and session are scalars indicating the subject and session of
%   the data files to be used. Filenames '3DS-<subject>-<session>.mat' and
%   '3DS-<subject>-<session>-eyepos.mat' are assumed in the current path.
%
% Outputs:
%   - h is a struct containing the handles for various items on the
% figure.  It can be used to get or change properties of the figures.  For 
% example, type 'get(h.Plot)' to get the properties of the main movie.
%
% Created 8/23/10 by DJ.
% Updated 8/25/10 by DJ - now only 1 click on bottom plot will select time
% Updated 9/1/10 by DJ - object visibility on time plot, target color
%    coding, x is now an input
% Updated 9/14/10 by DJ - added record_time to work with time offsets.
% Updated 5/16/11 by DJ - set any pupilsize NaN values to 1 to avoid error
% Updated 8/9/11 by DJ - added eye calibration
% Updated 2/8/13 by DJ - added saccade plotting, Jump-to-saccade buttons
% Updated 5/6/13 by DJ - added subject/session input format

% -------- INPUTS -------- %
if length(samples)==1 && length(pupilsize)==1    
    subject = samples;
    session = pupilsize;
    fprintf('Loading data from Subject %d, session %d...',subject,session);
    foo = load(sprintf('3DS-%d-%d-eyepos',subject,session));
    samples = foo.eyepos;
    pupilsize = foo.pupilsize;
    foo = load(sprintf('3DS-%d-%d',subject,session));
    x = foo.x;
    clear foo
    disp('Done!')
end

% -------- SETUP -------- %
fs = 1000; % sampling frequency of eyelink (Hz)
time = (1:length(samples))/fs; % time since beginning of eyelink sampling
ps_reg = 50/nanmax(pupilsize); % factor we use to regularize pupil size
screen_res = [1024, 768]; % size of screen in pixels (horizontal, vertical)
if nargin<4
    t_start = 1/fs;
end
% Parse inputs
object_limits = x.eyelink.object_limits;
record_time = x.eyelink.record_time; % time 'offset' - the time when eyelink started recording.
visible_times = GetVisibleTimes(x,'eeg');
isTarget = strcmp('TargetObject',{x.objects(:).tag});
colors = 'br';
% Get saccade info
saccade_start_pos = x.eyelink.saccade_start_positions;
saccade_end_pos = x.eyelink.saccade_positions;
saccade_times = ([x.eyelink.saccade_start_times, x.eyelink.saccade_times]-record_time)/1000; %[start end] in s
iSaccadeToObj = find(ismember(x.eyelink.saccade_times,x.eyelink.saccade_events(:,1)));

% Fix NaN problems
pupilsize(isnan(pupilsize)) = 1; % set non-existent pupilsize measures (i.e. during blinks) to tiny dot size
samples = ApplyEyeCalibration(samples,x); % apply eye calibration if you haven't already

% -------- INITIAL PLOTTING -------- %
disp('Setting up figure...');
figure; % make a new figure
MakeFigureTitle(x.EDF_filename,false);
[~,iTime] = min(abs(time-t_start)); % the global index of the current time point - this is used throughout all functions

% Main eye plot
h.Plot = axes('Units','normalized','Position',[0.13 0.3 0.775 0.65],'ydir','reverse'); % set position
hold on;
rectangle('Position',[0 0 screen_res]);
h.Dot = plot(samples(iTime,1),samples(iTime,2),'k.','MarkerSize',pupilsize(iTime)*ps_reg);
axis(h.Plot,[-200 screen_res(1)+200 -200 screen_res(2)+200]);
h.Rect = [];
h.iSac = [];
h.Saccade = [];

title(sprintf('t = %.3f s',time(iTime)));

% -------- TIME SELECTION PLOT SETUP -------- %
% 'Time plot' for selecting and observing current time
h.Time = axes('Units','normalized','Position',[0.13 0.1 0.775 0.1],'Yticklabel',''); % set position
hold on
% Plot horizontal lines when objects were visible
visible_objects = visible_times(:,1); % object numbers
for i=1:size(visible_times,1);
    obj = visible_objects(i); % the object that's visible
    plot((visible_times(i,2:3) - x.eeg.record_time)/x.eeg.eventsamplerate,...   % plot a horizontal line at the time this object was visible
        [obj obj],colors(isTarget(obj)+1),'ButtonDownFcn',@time_callback);      % color it as target/distractor and make it clickable
end
% Plot first saccade to each object
plot((x.eeg.saccade_events(:,1)-x.eeg.record_time)/x.eeg.eventsamplerate, ...
    x.eeg.saccade_events(:,2),'k+','ButtonDownFcn',@time_callback);
% Plot eye position
% normalized_samples = samples(:,1)/max(abs(samples(:,1)))*(max(visible_objects)-1);
% normalized_samples = normalized_samples - max(normalized_samples);
% plot(time,normalized_samples,'g','ButtonDownFcn',@time_callback); % plot eye position
% Annotate plot
plot([0 time(end)],[0 0],'k','ButtonDownFcn',@time_callback); % plot separation between plots
xlim(h.Time,[0 time(end)]);
xlabel('time (s)');
ylabel('eye pos   |   objects'); % Top section is object visibility, bottom section is eye x position
h.Line = plot(time([iTime iTime]), get(gca,'YLim'),'k','linewidth',2); % Line indicating current time
set(h.Time,'ButtonDownFcn',@time_callback) % this must be called after plotting, or it will be overwritten

% -------- GUI CONTROL SETUP -------- %
disp('Making GUI controls...')
h.Play = uicontrol('Style','togglebutton',...
                'String','Play',...
                'Units','normalized','Position',[.45 .2 .1 .05],...
                'Callback',@play_callback); % play button
h.Speed = uicontrol('Style','slider',...
                'Min',1,'Max',100,'Value',1,'SliderStep',[.05 .2],...
                'Units','normalized','Position',[.45 .25 .1 .025]); % speed slider
h.SacBack = uicontrol('Style','pushbutton',...
                'String','Sac <',...
                'Units','normalized','Position',[.25 .2 .1 .05],...
                'Callback',@sacback_callback); % back-to-last-saccade button
h.Back = uicontrol('Style','pushbutton',...
                'String','<',...
                'Units','normalized','Position',[.38 .2 .05 .05],...
                'Callback',@back_callback); % back button
h.Fwd = uicontrol('Style','pushbutton',...
                'String','>',...
                'Units','normalized','Position',[.57 .2 .05 .05],...
                'Callback',@fwd_callback); % forward button            
h.SacFwd = uicontrol('Style','pushbutton',...
                'String','> Sac',...
                'Units','normalized','Position',[.65 .2 .1 .05],...
                'Callback',@sacfwd_callback); % fwd-to-next-saccade button
              
disp('Done!')
    
% -------- SUBFUNCTIONS -------- %
function redraw() % Update the line and topoplot
    % Check that iTime is within allowable bounds
    if iTime<1 iTime=1;
    elseif iTime>numel(time) iTime = numel(time);
    end
    % Adjust plots
    set(h.Line,'XData',time([iTime iTime]));
    axes(h.Plot);
    if isnan(pupilsize(iTime))
        set(h.Dot,'MarkerSize',0.1)
    else
        set(h.Dot,'XData',samples(iTime,1),'YData',samples(iTime,2),'MarkerSize',pupilsize(iTime)*ps_reg);
    end
    
    % Plot Saccade
    saccade = find(time(iTime)>saccade_times(:,1) & time(iTime)<saccade_times(:,2),1);
    if ~isequal(saccade,h.iSac)
        delete(h.Saccade);
        h.Saccade = [];
        h.iSac = [];
    end
    if ~isempty(saccade) && isempty(h.Saccade)
        if any(iSaccadeToObj==saccade)
            saccolor = 'r.-';
        else
            saccolor = 'b.-';
        end            
        h.Saccade = plot([saccade_start_pos(saccade,1) saccade_end_pos(saccade,1)], [saccade_start_pos(saccade,2) saccade_end_pos(saccade,2)], saccolor);        
        h.iSac = saccade;    
    end    
    
    % Find object limits and plot them
    delete(h.Rect); % out with the old
    h.Rect = [];
    iLims = find(object_limits(:,1)>iTime+record_time-40 & object_limits(:,1)<iTime+record_time); % visible within the last 20 ms
    obj = unique(object_limits(iLims,2));
    for i=1:numel(obj)
        h.Rect(i) = rectangle('Position',object_limits(iLims(find(object_limits(iLims,2)==obj(i),1,'last')),3:6),'EdgeColor',colors(isTarget(obj(i))+1)); % in with the new
    end
    % Update title
    title(sprintf('t = %.3f s\n Object #%s, Saccade #%d',time(iTime),num2str(obj),saccade)); % use num2str in case there are multiple objects visible
end

function time_callback(hObject,eventdata) % First mouse click on the Time plot brings us here
    cp = get(h.Time,'CurrentPoint'); % get the point(s) (x,y) where the person just clicked
    x = cp(1,1); % choose the x value of one point (the x values should all be the same).
    iTime = find(time>=x,1); % find closest time to the click
    redraw; % update line and topoplot
end    

function play_callback(hObject,eventdata)
    % Get button value
    button_is_on = get(hObject,'Value') == get(hObject,'Max');
    % Keep incrementing and plotting
    while button_is_on && iTime < numel(time) %until we press pause or reach the end
        iTime=iTime+floor(get(h.Speed,'Value'));
        redraw;
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

function sacback_callback(hObject,eventdata)
    saccade = find(saccade_times(:,2)<time(iTime),1,'last');
    iTime = find(time>saccade_times(saccade,1),1);    
    redraw; % update line and topoplot
end

function sacfwd_callback(hObject,eventdata)
    saccade = find(saccade_times(:,1)>time(iTime),1,'first');
    iTime = find(time>saccade_times(saccade,1),1);    
    redraw; % update line and topoplot
end


end %function MakeEyeMovie



