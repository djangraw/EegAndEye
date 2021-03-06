function hOut = PlotEyeErps_MultiSession(subject,sessions,interpolateBlinks,epochRange)

% Combines the eye position data from multiple sessions and uses them as
% input to PlotEyeErps().
%
% hOut = PlotEyeErps_MultiSession(subject,sessions,interpolateBlinks)
%
% INPUTS:
% -subject is the subject number.
% -sessions is a vector of session numbers.
% The filename '3DS-<subject>-<sessions(i)>' should be a 3DS data file
% imported using Import_3DS_Data_v3. '3DS-<subject>-<sessions(i)>-eyepos'
% should be the raw eye position data imported using find_events().
% -interpolateBlinks is a binary value indicating whether we should
% interpolate data during blinks. [default: true]
%
% OUTPUTS:
% -hOut is a struct containing handles of various parts of the plot (see
% PlotEyeErps) as well as hOut.data, a struct matching the output of
% GetEyeErps but with all the indicated sessions included, and 
% hOut.FigureTitle, the text at the top of the figure.
%
% Created 3/29/11 by DJ.
% Updated 8/9/11 by DJ - added interpolateBlinks input
% Updated 11/3/11 by DJ - added epochRange input to avoid end-of-session 
%   error, and swithed order of ApplyEyeCalibration to avoid blink warning

if nargin<3 || isempty(interpolateBlinks)
    interpolateBlinks = true;
end
if nargin<4
    epochRange = [-500 3000];
end

fprintf('---Creating eye ERPs for Subject %d, %d sessions\n',subject,numel(sessions));
BigData = struct('epochSamples',[],'epochTimes',[],'targetEventSessions',[],...
    'targetEventTimes',[],'targetEyeEpochs',[],'targetPupEpochs',[],...
    'targetObjEpochs',[],'distractorEventSessions',[],...
    'distractorEventTimes',[],'distractorEyeEpochs',[],...
    'distractorPupEpochs',[],'distractorObjEpochs',[],'screen_res',[]);
for i=1:numel(sessions)
    fprintf('---loading file %d of %d...\n',i,numel(sessions));
    % load data
    filename = sprintf('3DS-%d-%d',subject,sessions(i));
    load(filename) % get x
    load([filename '-eyepos']); % get eyepos, pupilsize
    
    % calibrate eye position
    eyepos_calibrated = eyepos;%ApplyEyeCalibration(eyepos,x);
    
    if interpolateBlinks
        % interpolate around blinks
        eyepos = InterpolateBlinks(eyepos,x.eyelink.record_time-1+(1:length(eyepos)),x);
        pupilsize = InterpolateBlinks(pupilsize,x.eyelink.record_time-1+(1:length(pupilsize)),x);
    end
    
    % Get the ERPs!
    data = GetEyeErps(eyepos_calibrated,pupilsize,x,epochRange);
    
    % add to BigData
    BigData.epochSamples = data.epochSamples;
    BigData.epochTimes = data.epochTimes;
    BigData.targetEventSessions = [BigData.targetEventSessions; repmat(sessions(i),size(data.targetEventTimes))];
    BigData.targetEventTimes = [BigData.targetEventTimes; data.targetEventTimes];
    BigData.targetEyeEpochs = [BigData.targetEyeEpochs; data.targetEyeEpochs];
    BigData.targetPupEpochs = [BigData.targetPupEpochs; data.targetPupEpochs];
    BigData.targetObjEpochs = [BigData.targetObjEpochs; data.targetObjEpochs];
    BigData.distractorEventSessions = [BigData.distractorEventSessions; repmat(sessions(i),size(data.distractorEventTimes))];
    BigData.distractorEventTimes = [BigData.distractorEventTimes; data.distractorEventTimes];
    BigData.distractorEyeEpochs = [BigData.distractorEyeEpochs; data.distractorEyeEpochs];
    BigData.distractorPupEpochs = [BigData.distractorPupEpochs; data.distractorPupEpochs];
    BigData.distractorObjEpochs = [BigData.distractorObjEpochs; data.distractorObjEpochs];
    BigData.screen_res = data.screen_res;
end
    
disp('---Plotting eye ERPs...')
% Plot eye erps from BigData
hOut = PlotEyeErps(BigData);
hOut.FigureTitle = MakeFigureTitle(sprintf('Subject %d Eye Position',subject));

disp('---Done!')