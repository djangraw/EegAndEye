function TestFlightSimFile(subject,session)
%
% Created 9/4/14 by DJ.

% Declare parameters
% experiment = 'FlightSim';
experiment = 'FlightSim_1pt0';
% subject = 1;
% session = 9;

% Import data
NEDE_ImportData(experiment,subject,session);
% Load resulting file
filename = sprintf('%s-%d-%d.mat',experiment,subject,session);
load(filename);

% Import eye samples
cd('samples');
GetEyeSamples(subject,session,experiment,0);
load(sprintf('%s-%d-%d-eyepos.mat',experiment,subject,session));
cd ..


%% Extract position
ringpos = nan(numel(x.objects),2);
for i=1:numel(x.objects)
    ringpos(i,:) = [x.objects(i).position(:,2), x.objects(i).elevation]; % [z, y]
end
% Extract size
ringsize = nan(numel(x.objects),2);
for i=1:numel(x.objects)
    ringsize(i,:) = x.objects(i).rotation([3 2]); % [z, y]
end

% extract subject path
subjpos = [x.events.camera.position(:,2), x.events.camera.elevation]; % z,y
% interpolate to get constant sampling frequency
zpos = linspace(min(subjpos(:,1)),max(subjpos(:,1)),size(subjpos,1));
ypos = interp1(subjpos(:,1),subjpos(:,2),zpos);
% replace subject path with interpolated version
subjpos = [zpos',ypos'];
isInFlight = subjpos(:,1)>0; % post-takeoff, when subject had control

% extract drift times & positions
iDriftMsg = find(strncmp('Drift',x.events.message.text,5));
[drift,tDrift] = deal(zeros(1,numel(iDriftMsg)));
driftPos = zeros(numel(iDriftMsg),2);
for i=1:numel(iDriftMsg)
    drift(i) = str2double(x.events.message.text{iDriftMsg(i)}(8:end));
    tDrift(i) = x.events.message.time(iDriftMsg(i));
    [~,iDriftPos] = min(abs(x.events.camera.time-tDrift(i)));
    driftPos(i,:) = [x.events.camera.position(iDriftPos,2),x.events.camera.elevation(iDriftPos)];
end

% extract controls
iCtrlMsg = find(strncmp('Control',x.events.message.text,5));
[ctrl,tCtrl] = deal(zeros(1,numel(iCtrlMsg)));
ctrlPos = zeros(numel(iCtrlMsg),2);
iParen = find(x.events.message.text{iCtrlMsg(1)}=='(');
for i=1:numel(iCtrlMsg)
    iComma = find(x.events.message.text{iCtrlMsg(i)}==',');
    ctrl(i) = str2double(x.events.message.text{iCtrlMsg(i)}(iParen+1:iComma-1));
    tCtrl(i) = x.events.message.time(iCtrlMsg(i));
    [~,iCtrlPos] = min(abs(x.events.camera.time-tCtrl(i)));
    ctrlPos(i,:) = [x.events.camera.position(iCtrlPos,2),x.events.camera.elevation(iCtrlPos)];
end

% Extract pupil size samples info
startTime = x.events.recording.time_start;
eyeTimes = (1:length(eyepos))+startTime-1;

blinkExpand = [-50 50];
blinkRanges = [x.events.blink.time_start+blinkExpand(1),x.events.blink.time_end+blinkExpand(2)];
pup_interp = InterpolateBlinks(pupilsize,eyeTimes,blinkRanges);

eye_z = interp1(x.events.camera.time,x.events.camera.position(:,2),eyeTimes);


%% Plot
% --- Set up plot
clf; 
subplot(4,1,1); hold on;
ylim([-4 4]);
% Plot subject path
% plot(ctrlPos(:,1),ctrl);
plot(ctrlPos(:,1),SmoothData(ctrl,5,'half'));
% Plot ideal path (through ring centers)
plot(ringpos(:,1),zeros(size(ringpos,1),1),'r.-');
% plot ground
plot(get(gca,'xlim'),[0 0],'k--');
% plot drift change positions
PlotVerticalLines(driftPos(:,1),'g--');
% plot drift angles
for i=1:size(driftPos,1)
    if i<size(driftPos,1)
        plot([driftPos(i,1),driftPos(i+1,1)], [drift(i) drift(i)],'g');
    else
        plot([driftPos(i,1),subjpos(end,1)], [drift(i) drift(i)],'g');
    end
%     plot(driftPos(i,1)+[0 200], [0 200*tan(drift(i)*pi/180)],'g-');
end


% Annotate plot
xlim([min(subjpos(:,1)), max(subjpos(:,1))]);
xlabel('Z position (m)');
ylabel('Joystick Input (A.U.)');
title(show_symbols(filename))
legend('Joystick Input','Ring Positions','Zero','Drift Change','Drift Amplitude','Location','EastOutside');





% --- Plot subject path
subplot(4,1,2); hold on;
plot(subjpos(~isInFlight,1),subjpos(~isInFlight,2),'b--');
plot(subjpos(isInFlight,1),subjpos(isInFlight,2));
% Plot ideal path (through ring centers)
plot(ringpos(:,1),ringpos(:,2),'r.-');
% Plot top edge
plot(ringpos(:,1),ringpos(:,2)+ringsize(:,2)/2,'r:');
% plot ground
plot(get(gca,'xlim'),[0 0],'k');
% plot drift change positions
PlotVerticalLines(driftPos(:,1),'g--');
% % plot drift angles
% for i=1:size(driftPos,1)
%     plot(driftPos(i,1)+[0 200], [0 200*tan(drift(i)*pi/180)],'g-');
% end

% plot bottom edge (here to fix legend)
plot(ringpos(:,1),ringpos(:,2)-ringsize(:,2)/2 ,'r:');
% Plot rings
for i=1:numel(x.objects)
    plot([ringpos(i,1), ringpos(i,1)], [ringpos(i,2)-ringsize(i,2)/2, ringpos(i,2)+ringsize(i,2)/2],'r-');        
end

% Annotate plot
xlim([min(subjpos(:,1)), max(subjpos(:,1))]);
xlabel('Z position (m)');
ylabel('Elevation (m)');
title(show_symbols(filename))
legend('Subject Path: takeoff','Subject Path: flight','Rings','Boundaries','Ground','Drift Change','Location','EastOutside');

% --- Set up PIO plot
subplot(4,1,3); hold on;
idealPath = [subjpos(isInFlight,1),interp1(ringpos(:,1),ringpos(:,2),subjpos(isInFlight,1))];
pathError = [zeros(sum(~isInFlight),1); -(idealPath(:,2)-subjpos(isInFlight,2))];

hilbPathError = hilbert(pathError);

% plot(subjpos(:,1),abs(pathError));
plot(subjpos(~isInFlight,1),abs(hilbPathError(~isInFlight)),'b--');
plot(subjpos(isInFlight,1),abs(hilbPathError(isInFlight)));
plot(ringpos(:,1),zeros(size(ringpos,1),1),'r.-');
plot(ringpos(:,1),ringsize(:,2)/2,'r:');
% plot drift change positions
PlotVerticalLines(driftPos(:,1),'g--');
% % plot drift angles
% for i=1:size(driftPos,1)
%     plot(driftPos(i,1)+[0 200], [0 200*tan(drift(i)*pi/180)],'g-');
% end
    
% plot individual rings
for i=1:numel(x.objects)
    plot([ringpos(i,1), ringpos(i,1)], [0 ringsize(i,2)/2],'r-');        
end


xlim([min(subjpos(:,1)), max(subjpos(:,1))]);
xlabel('Z position (m)');
ylabel('PIO Amplitude');
% ylabel('Hilbert Transform of error (m)');
legend('hilb(error=0): takeoff','hilb(error): flight','Rings','Boundary','Drift Change','Location','EastOutside');


% --- Set up pupil size plot
subplot(4,1,4); hold on;
% Plot subject path
% plot(ctrlPos(:,1),ctrl);
plot(eye_z,pup_interp-mean(pup_interp));
% Plot ideal path (through ring centers)
plot(ringpos(:,1),zeros(size(ringpos,1),1),'r.-');
% plot ground
plot(get(gca,'xlim'),[0 0],'k--');
% plot drift change positions
PlotVerticalLines(driftPos(:,1),'g--');
% % plot drift angles
% for i=1:size(driftPos,1)
%     plot(driftPos(i,1)+[0 200], [0 200*tan(drift(i)*pi/180)],'g-');
% end
% Annotate plot
xlim([min(subjpos(:,1)), max(subjpos(:,1))]);
xlabel('Z position (m)');
ylabel('Pupil Size (A.U.)');
% title(show_symbols(filename))
legend('Pupil Size','Ring Positions','Zero','Drift Change','Location','EastOutside');





% Standardize plot widths
plotpos = zeros(4,4);
for i=1:4
    plotpos(i,:) = get(subplot(4,1,i),'Position');
end
minwidth = min(plotpos(:,3));
plotpos(:,3) = minwidth;
for i=1:4    
    set(subplot(4,1,i),'Position',plotpos(i,:));    
end
