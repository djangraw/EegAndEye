function PlotObjectLifetimes(x)

% PlotObjectLifetimes(x)
%
% - Plots the lifetime of each object as a line, w/ visible times in bold,
% saccades to that object as +'s, and button presses as vertical lines.
% - INPUT x is ImportData's x, a struct containing info about the session 
% and the objects in it.
%
% Created 6/24/10 by DJ.
% Updated 7/26/10 by DJ - new input format, added button times
% Updated 7/29/10 by DJ - changed events field back to eyelink
% Updated 8/19/10 by DJ - introduced saccade_events as a separate field
% Updated 12/5/13 by DJ - include leader events

% SETUP
GetNumbers;
nObjects = numel(x.objects);
Ts = x.eeg.start_time; % time when the recording started
Te = x.eeg.end_time; % time when the trial ended
lifetimes = x.eeg.object_lifetimes; % create and destroy times of each object
vision_events = x.eeg.object_events; % when objects entered, exited the scene
saccade_events = x.eeg.saccade_events; % usually the first saccade to each object
button_times = x.eeg.button_times; % when a target button was pressed
brake_button_times = x.eeg.brake_button_times; % when a brake button was pressed
leader_slow_times = x.eeg.leader_slow_times; % when the leading car slowed down
fs = x.eeg.eventsamplerate; % sampling rate of eeg
visible_times = GetVisibleTimes(x);

% set up axis
cla;
hold on
colors = 'brkmg'; % colors for distractor, target, saccade, target button and brake button
markers = {'b:','r:','r','k+','m-','g--','g-'};

for i=1:nObjects
    % Get object info
    object = x.objects(i);
    isTarget = strcmp(object.tag,'TargetObject');
    
    % Plot lifetime as a thin line
    plot((lifetimes(i,:)-Ts)/fs, [i i],markers{isTarget+1});
    
    VT_isthisobj = visible_times(:,1)==i;
    enter_exit = visible_times(VT_isthisobj,2:3);
    for j=1:size(enter_exit,1)
        plot((enter_exit(j,:) - Ts)/fs,[i,i],markers{isTarget+1}(1),'linewidth',2);
    end
    
    % Plot saccade times as black +'s
    if ~isempty(saccade_events)
        x_times = saccade_events(saccade_events(:,2)==(Numbers.SACCADE_TO+i),1);
        plot((x_times - Ts)/fs,repmat(i,1,numel(x_times)),markers{4},'MarkerSize',12);
    end
end

% Plot button and brake times as vertical dotted lines
for i=1:numel(button_times)
    plot(([button_times(i) button_times(i)]-Ts)/fs, [0 nObjects+1], markers{5});
end 

for i=1:numel(leader_slow_times)
    plot(([leader_slow_times(i) leader_slow_times(i)]-Ts)/fs, [0 nObjects+1], markers{6});
end 

for i=1:numel(brake_button_times)
    plot(([brake_button_times(i) brake_button_times(i)]-Ts)/fs, [0 nObjects+1], markers{7});
end 


% Annotate plot
if isfield(x,'EDF_filename')
    title(show_symbols(sprintf('Visibility of objects in %s',x.EDF_filename)));
else
    title('Visibility of objects in session');
end
xlim([0 (Te - Ts)/fs]);
ylim([0 nObjects+1]);
xlabel('time (s)');
ylabel('object number');
MakeLegend(markers,{'Distractor','Target','Visible','1st Saccade To Obj','TargButton','LeaderSlow','BrakeButton'},[1 1 2 1 1 1,1]);