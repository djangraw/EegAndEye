function hOut = PlotEyeAndEeg(data,trials,EEG,isTargetTrial,offsets)

% Makes a figure/UI for scrolling through eye position and EEG data 
% side-by-side for a specified group of trials.
%
% hOut = PlotEyeAndEeg(data,trials,EEG,offsets)
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
%   - trials is a vector of trial numbers in data.
%   - EEG is the corresponding EEGLAB dataset.
%   - isTargetTrial is a binary value that indicates whether the trials
%   are target trials. [default: infer from EEG setname]
%   - offsets (not yet functional) is a vector of time offsets (in ms) for each
%   trial to line them up in a specific way.
%
% Outputs:
%   - hOut is a struct containing the handles for various items on the
% figure.  It can be used to get or change properties of the figures.  For 
% example, type 'get(h.EyePlot)' to get the properties of the eye movie.
%
% Created 5/12/11 by DJ based on PlotEyeErps.m.
% Updated 5/13/11 by DJ - added default for isTargetTrial
% Updated 7/21/11 by DJ - added pupil size plot, separate trial saccades

% Read from EEG setname whether this is a target EEG file
if nargin<4
    if isempty(strfind(EEG.setname,'targ'))
        isTargetTrial = 0;
    else
        isTargetTrial = 1;
    end
end

% -------- UNPACK -------- %
epochSamples = data.epochSamples;
screen_res = data.screen_res;

% -------- CROP EYE DATA -------- %
if isTargetTrial
    data.targetEventSessions = data.targetEventSessions(trials);
    data.targetEventTimes = data.targetEventTimes(trials);
    data.targetEyeEpochs = data.targetEyeEpochs(trials,:);
    data.targetPupEpochs = data.targetPupEpochs(trials,:);
    data.targetObjEpochs = data.targetObjEpochs(trials,:);
else
    data.targetEventSessions = data.distractorEventSessions(trials);
    data.targetEventTimes = data.distractorEventTimes(trials);
    data.targetEyeEpochs = data.distractorEyeEpochs(trials,:);
    data.targetPupEpochs = data.distractorPupEpochs(trials,:);
    data.targetObjEpochs = data.distractorObjEpochs(trials,:);
end


% -------- EXTRACT EEG DATA -------- %
trialTimes = data.targetEventTimes;
eegTrials = TimeToEpochNumber(EEG,trialTimes);
eegData = mean(EEG.data(:,:,eegTrials),3);
eegTimes = EEG.times;
eegLocs = EEG.chanlocs;
pupilTimes = data.epochTimes;
pupilData = data.targetPupEpochs;
% find color limits
color_max = max(eegData(:));
colorLimits = [-color_max color_max]; %topoplot color limits (in uV for an ERP, ERP ylabel displays these units)
colorLimits = [-10 10];

% -------- INITIAL PLOTTING -------- %
disp('Setting up figure...');
figure; % make a new figure
% global iStart iEnd; % the global index of the current time point - this is used throughout all functions
iStart = 1; 
iEnd = 2;

% Eye plot
h.EyePlot = axes('Units','normalized','Position',[0.13 0.3 0.35 0.65]); % set position
rectangle('Position',[0 0 screen_res]);
title('Eye Position');
axis(h.EyePlot,[-200 screen_res(1)+200 -200 screen_res(2)+200]);

% EEG plot
h.EegPlot = axes('Units','normalized','Position',[0.53 0.3 0.35 0.65]); % set position
zeroTime = eegTimes==0;
topoplot(eegData(:,zeroTime),eegLocs,'MapLimits',colorLimits,'conv','on'); % initial topoplot
colorbar; % add colorbar (this only needs to be called once)
title('Mean Voltage');

% -------- TIME SELECTION PLOT SETUP -------- %
% 'Time plot' for selecting and observing current time
h.Time = axes('Units','normalized','Position',[0.13 0.1 0.775 0.1],'Yticklabel','','ButtonDownFcn',@time_callback); % set position
% set(h.Time,'Xticklabel',get(gca,'XTick'));
hold on;

lineSize = 1000;
% Plot eeg
% plot(eegTimes,eegData);
% plot pupil size
for i=1:size(pupilData,1)
    plot(pupilTimes,pupilData(i,:) - pupilData(i,pupilTimes==0)+ lineSize*(i-0.5));
end

% Saccade times
saccade_times = GetEpochSaccades(EEG);
saccade_times = saccade_times(eegTrials);
for i=1:numel(saccade_times)
    for j = 1:numel(saccade_times{i})
        plot(repmat(saccade_times{i}(j),1,2), [0 lineSize] + lineSize*(i-1),'k');
    end
end

% Annotate plot
% xlim([epochSamples(1) epochSamples(end)]);
xlim([-500 2000]);
ylim([0 lineSize*size(pupilData,1)])
xlabel('time (samples)');
plot([0 0],get(gca,'YLim'),'k');

% Plot time selection rectangle
h.Rect = imrect(h.Time,[iStart, min(get(h.Time,'YLim')),iEnd-iStart,range(get(h.Time,'YLim'))],...
    'PositionConstraintFcn',@roi_constraint);
addNewPositionCallback(h.Rect,@roi_callback); % Line indicating current time
set(h.Time,'ButtonDownFcn',@time_callback); % this must be called after plotting, or it will be overwritten
title('Drag rectangle to set limits of plot')      
     


% -------- SEND INFO TO FIGURE USERDATA -------- %
data.eegData = eegData;
data.eegTimes = eegTimes;
data.eegLocs = eegLocs;
data.colorLimits = colorLimits;
h.data = data;
set(gcf,'UserData',h);
hOut = h;
disp('Done!')
end %function 
    
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
    eegData = h.data.eegData;
    eegTimes = h.data.eegTimes;
    eegLocs = h.data.eegLocs;

    
    % Adjust plots
    title(h.Time,sprintf('From t=%.1f to t=%.1f ms',h.data.epochTimes(epochSamples==iStart),h.data.epochTimes(epochSamples==iEnd)));

    
    % Plot lines for target and distractor trials
    axes(h.EyePlot); cla; hold on;   
    for i=1:size(tEye,1) % for each target trial
        % plot object position
        trialObjEnd = tObj{i,end};
        if ~isempty(trialObjEnd)
            rectangle('Position',[trialObjEnd],'EdgeColor','k');
        end
        % plot eye position
        trialEye = [tEye{i,:}];
        plot(trialEye(1,:),trialEye(2,:),'b-','UserData',...
            sprintf('clicked target trial %d (session %d, time %.3f)',i,tSess(i),tTime(i)),...
            'ButtonDownFcn','disp(get(gco,''UserData''))');
        plot(trialEye(1,end),trialEye(2,end),'bo','UserData',...
            sprintf('clicked target trial %d (session %d, time %.3f)',i,tSess(i),tTime(i)),...
            'ButtonDownFcn','disp(get(gco,''UserData''))');
    end
    rectangle('Position',[0 0 screen_res]);
   
    % Plog EEG
    axes(h.EegPlot);
    eegSamples = eegTimes>=h.data.epochTimes(iStart-epochSamples(1)) & eegTimes <= h.data.epochTimes(iEnd-epochSamples(1));
    topoplot(mean(eegData(:,eegSamples),2),eegLocs,'MapLimits',h.data.colorLimits,'conv','on');
    
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

end    






