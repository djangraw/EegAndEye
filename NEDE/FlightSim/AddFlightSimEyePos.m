function AddFlightSimEyePos(filename,x,pupilsize)

% Created 9/23/14 by DJ.

% set up
blinkExpand = [-50 50];
blinkRanges = [x.events.blink.time_start+blinkExpand(1),x.events.blink.time_end+blinkExpand(2)];
nTrials = length(x.events.trial.time_start);

% Create
[pup_cell,pupPos_cell] = deal(cell(1,nTrials));
for iTrial = 1:nTrials
    if mod(iTrial/nTrials*100,10)==0
        fprintf('Trial %d/%d...\n',iTrial,nTrials);
    end
    tStart = x.events.trial.time_start(iTrial);
    tEnd = x.events.trial.time_end(iTrial);
    
    startTime = x.events.recording.time_start(iTrial);         
    tEye = (1:length(pupilsize{iTrial}))+startTime-1;
    pup_interp = InterpolateBlinks(pupilsize{iTrial},tEye,blinkRanges);
    eye_z = interp1(x.events.camera.time,x.events.camera.position(:,2),tEye);
    
    pupPos_cell{iTrial} = eye_z(tEye>=tStart & tEye<=tEnd);
    pup_cell{iTrial} = pup_interp(tEye>=tStart & tEye<=tEnd);
end


fprintf('loading...\n');
load(filename);

fprintf('Saving...\n');
save(filename,'subjpos_cell','isInFlight_cell','ctrlPos_cell','ctrl_cell','pupPos_cell','pup_cell','zlimits','ringpos','ringsize','filename','nTrials');
fprintf('Done!\n');