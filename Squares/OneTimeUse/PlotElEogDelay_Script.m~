%% S6, N2: Load data
load sq-6-2
load sq-6-2-eyepos
eegdata = readEEG_b4preprocessing('sq-6.001.b',1,[],0); % ch. 1 = heogr
eegevents = get_sensorium_events('sq-6.001.b'); % sync events

% Set times
t = 8900:10900; % eeg times of interest
delay = round(mean(x.sync.eyelink(1:5)-eegevents(1:5,1)));
starttime = 86242550; % from .asc events file

%% S6, N7: Load data
load sq-6-7
load sq-6-7-eyepos
eegdata = readEEG_b4preprocessing('sq-6.006.b',1,[],0); % ch. 1 = heogr
eegevents = get_sensorium_events('sq-6.006.b'); % sync events

% Set times
t = 6100:8100; % eeg times of interest
delay = round(mean(x.sync.eyelink(1:5)-eegevents(1:5,1)));
starttime = 88616557; % from .asc events file

%% S7, N2: Load data
load sq-7-2
load sq-7-2-eyepos
eegdata = readEEG_b4preprocessing('sq-7.001.b',1,[],0); % ch. 1 = heogr
eegevents = get_sensorium_events('sq-7.001.b'); % sync events

% Set times
t = 12400:14400; % eeg times of interest
delay = round(mean(x.sync.eyelink(1:5)-eegevents(1:5,1)));
starttime = 29478; % from .asc events file

%% S7, N8: Load data
load sq-7-8
load sq-7-8-eyepos
eegdata = readEEG_b4preprocessing('sq-7.007.b',1,[],0); % ch. 1 = heogr
eegevents = get_sensorium_events('sq-7.007.b'); % sync events

% Set times
t = 9500:11500; % eeg times of interest
delay = round(mean(x.sync.eyelink(1:5)-eegevents(1:5,1)));
starttime = 22573; % from .asc events file


%% S8, N2: Load data
load sq-8-2
load sq-8-2-eyepos
eegdata = readEEG_b4preprocessing('sq-8.001.b',1,[],0); % ch. 1 = heogr
eegevents = get_sensorium_events('sq-8.001.b'); % sync events

% Set times
t = 300:2300; % eeg times of interest
delay = round(mean(x.sync.eyelink(1:5)-eegevents(1:5,1)));
starttime = 26559; % from .asc events file

%% S8, N9: Load data
load sq-8-9
load sq-8-9-eyepos
eegdata = readEEG_b4preprocessing('sq-8.008.b',1,[],0); % ch. 1 = heogr
eegevents = get_sensorium_events('sq-8.008.b'); % sync events

% Set times
t = 4200:6200; % eeg times of interest
delay = round(mean(x.sync.eyelink(1:5)-eegevents(1:5,1)));
starttime = 36335; % from .asc events file

%% Do it up

% Plot
clf
cla
hold on
plot(t/1000,eyepos(t+delay-starttime,1)')
plot(t/1000,eegdata(1,t)'*1e6,'g')
PlotVerticalLines(eegevents(1:3)/1000,'r')
oksaccades = find(x.saccade.start_time-delay>t(1) & x.saccade.start_time-delay<t(end));
PlotVerticalLines((x.saccade.start_time(oksaccades)-delay)/1000,'k--')

% Annotate plot
xlim([t(1) t(end)]/1000)
xlabel('eeg clock time (s)')
ylabel('voltage/position (a.u.)')
title(sprintf('S%d, Session %d Delay',x.subject,x.session))
MakeLegend({'b' 'g' 'r' 'k--'},{'eye x position','HEOGR','sync event','saccade start time'});



















