function filename = GetFlightSimBehavior(x,pupilsize)

% GetFlightSimBehavior(x,pupilsize)
% 
% INPUTS:
%
% Created 9/22/14 by DJ based on TestFlightSimFile.m.


%% Extract position

fprintf('%s: Extracting Ring Info...\n',datestr(now,16))
ringpos = nan(numel(x.objects),2);
for i=1:numel(x.objects)
    ringpos(i,:) = [x.objects(i).position(:,2), x.objects(i).elevation]; % [z, y]
end
% Extract size
ringsize = nan(numel(x.objects),2);
for i=1:numel(x.objects)
    ringsize(i,:) = x.objects(i).rotation([3 2]); % [z, y]
end

fprintf('%s: Extracting Path Info...\n',datestr(now,16))
% extract subject path
subjpos = [x.events.camera.position(:,2), x.events.camera.elevation]; % z,y
isInFlight = subjpos(:,1)>0; % post-takeoff, when subject had control
tSubj = x.events.camera.time; % not quite right, but ok.

fprintf('%s: Extracting Controls Info...\n',datestr(now,16))
% extract controls
iCtrlMsg = find(strncmp('Control',x.events.message.text,5));
[ctrl,tCtrl] = deal(zeros(1,numel(iCtrlMsg)));
ctrlPos = zeros(numel(iCtrlMsg),2);
iParen = find(x.events.message.text{iCtrlMsg(1)}=='(',1);
for i=1:numel(iCtrlMsg)
    iComma = find(x.events.message.text{iCtrlMsg(i)}==',',1);
    ctrl(i) = str2double(x.events.message.text{iCtrlMsg(i)}(iParen+1:iComma-1));
    tCtrl(i) = x.events.message.time(iCtrlMsg(i));
    [~,iCtrlPos] = min(abs(x.events.camera.time-tCtrl(i)));
    ctrlPos(i,:) = [x.events.camera.position(iCtrlPos,2),x.events.camera.elevation(iCtrlPos)];
end

% Interpolate blinks
if exist('pupilsize','var')
    blinkExpand = [-50 50];
    blinkRanges = [x.events.blink.time_start+blinkExpand(1),x.events.blink.time_end+blinkExpand(2)];
end


%%
nTrials = numel(x.events.trial.time_start);
[subjpos_cell,isInFlight_cell,ctrlPos_cell,ctrl_cell,pupPos_cell,...
    pup_cell,] = deal(cell(1,nTrials));
fprintf('%s: Parsing %d trials...\n',datestr(now,16),nTrials)
for iTrial=1:nTrials
    if mod(round(iTrial/nTrials*100),10) == 0
        fprintf('%d%% done...\n',round(iTrial/nTrials*100));
    end
    tStart = x.events.trial.time_start(iTrial);
    tEnd = x.events.trial.time_end(iTrial);    

    % get take out events for this trial

    % interpolate to get constant sampling frequency
    subjpos0 = subjpos(tSubj>=tStart & tSubj<=tEnd,:);
    zpos = linspace(min(subjpos0(:,1)),max(subjpos0(:,1)),size(subjpos0,1));
    ypos = interp1(subjpos0(:,1),subjpos0(:,2),zpos);
    % replace subject path with interpolated version
    subjpos1 = [zpos',ypos'];


    subjpos_cell{iTrial} = subjpos1;%subjpos(tSubj<=tStart & tSubj<=tEnd,:);
    isInFlight_cell{iTrial} = isInFlight(tSubj>=tStart & tSubj<=tEnd,:);
    ctrlPos_cell{iTrial} = ctrlPos(tCtrl>=tStart & tCtrl<=tEnd,:);
    ctrl_cell{iTrial} = ctrl(tCtrl>=tStart & tCtrl<=tEnd);

    % Extract pupil size samples info
    if exist('eyepos','var')       
        startTime = x.events.recording.time_start(iTrial);         
        tEye = (1:length(pupilsize{iTrial}))+startTime-1;
        pup_interp = InterpolateBlinks(pupilsize{iTrial},tEye,blinkRanges);
        eye_z = interp1(x.events.camera.time,x.events.camera.position(:,2),tEye);

        pupPos_cell{iTrial} = eye_z(tEye>=tStart & tEye<=tEnd);
        pup_cell{iTrial} = pup_interp(tEye>=tStart & tEye<=tEnd);
    end
end

zlimits = [min(subjpos(:,1)),max(subjpos(:,1))]; % for 

%% Save results

iHyphen = find(x.params.EDF_filename=='-',1);
filename = sprintf('%s-%d-all-Behavior.mat',x.params.EDF_filename(1:iHyphen-1),x.params.subject);
fprintf('%s: Saving as %s...\n',datestr(now,16),filename)
save(filename,'subjpos_cell','isInFlight_cell','ctrlPos_cell','ctrl_cell','pup_cell','pupPos_cell','zlimits','ringpos','ringsize','filename','nTrials');