function R = ClassifyWithEegAndDwellTime(ALLEEG,y,cvmode,offset,dwell,useEEG)

% Created 5/22/13 by DJ.
% Updated 7/11/13 by DJ - added support for as-is dwell feature, useEEG input
% Updated 8/26/13 by DJ - added y_EEG field to Results struct

if ~exist('offset','var')
    offset = 0;
end

if ~exist('dwell','var')
    %% Get dwell time
    % Calculate dwell tme
    disp('Getting dwell time...')
    [dwell_cell isToTarget_cell] = GetDwellTimes(y);
    dwell_vec = [dwell_cell{:}];
    isTarget_vec = [isToTarget_cell{:}];
    % crop to real trials
    dwell = dwell_vec(~isnan(dwell_vec));
    isTarget = isTarget_vec(~isnan(dwell_vec));

    %% Crop dwell trials
    % Reject same trials in dwell that you rejected in EEG
    if isfield(ALLEEG(1).etc,'rejectepoch')
        if isequal(ALLEEG(1).etc.rejectepoch,ALLEEG(2).etc.rejectepoch)
            fprintf('removing %d trials to match EEG record...\n',sum(ALLEEG(1).etc.rejectepoch));
            dwell(ALLEEG(1).etc.rejectepoch) = []; % reject both at once so trials still match up
            isTarget(ALLEEG(1).etc.rejectepoch) = [];
        else
            distractors = find(~isTarget);
            targets = find(isTarget);
            isBadDistractor = ALLEEG(1).etc.rejectepoch;
            isBadTarget = ALLEEG(2).etc.rejectepoch;
            fprintf('removing %d trials to match EEG record...\n',sum(isBadDistractor)+sum(isBadTarget));
            dwell([distractors(isBadDistractor), targets(isBadTargets)]) = []; % reject both at once so trials still match up
            isTarget([distractors(isBadDistractor), targets(isBadTargets)]) = [];
        end
    end
end

if ~exist('useEEG','var')
    useEEG = true;
end

%% Sort EEG trials
% put ALLEEG data into single matrix
data_0 = cat(3,ALLEEG(1).icaact,ALLEEG(2).icaact);
fwdModelData_0 = cat(3,ALLEEG(1).data,ALLEEG(2).data);
truth_0 = [zeros(1,ALLEEG(1).trials), ones(1,ALLEEG(2).trials)];
% find original trial order
epochID1 = zeros(1,ALLEEG(1).trials);
epochID2 = zeros(1,ALLEEG(2).trials);
for j=1:ALLEEG(1).trials
    epochID1(j) = ALLEEG(1).epoch(j).eventurevent{1}; 
end
for j=1:ALLEEG(2).trials
    epochID2(j) = ALLEEG(2).epoch(j).eventurevent{1}; 
end
% sort by trial order
[~,order] = sort([epochID1, epochID2],'ascend');
data = data_0(:,:,order);
fwdModelData = fwdModelData_0(:,:,order);
truth = truth_0(order);

% if ~isequal(truth,isTarget)
%     error('EEG and dwell-time truth vectors do not match!')
% end


%% Classify
% Declare time windows for RSVP classifier
if useEEG
    tWindowStart = 100:100:900;
else
    tWindowStart = [];
end
tWindowLength = 100;
% Convert to samples
timesWithOffset = ALLEEG(1).times - offset;
samplesWindowStart = round(interp1(timesWithOffset, 1:ALLEEG(1).pnts, tWindowStart));
samplesWindowLength = round(tWindowLength/1000*ALLEEG(1).srate);

% Run RSVP classifier
if size(dwell,1)==1
    dwellFeature = dwell'-mean(dwell); % transpose (make rows=trials) and de-mean
else
    dwellFeature = dwell;
end
% --- TEMP - USE JUST TOP 10 ICS AND NO DWELL TIME!
% dwellFeature = [];
% data = data(1:10,:,:);
% --- END TEMP
[y, w, v, fwdModel, y_EEG] = Run2LevelClassifier_nested(data,truth,samplesWindowLength,samplesWindowStart,cvmode,dwellFeature,fwdModelData);
% [y2, w2, v2, fwdModel2] = ...
% run_rsvp_classifier_rawdata(data,truth,samplesWindowLength,samplesWindowStart,cvmode,dwellFeature,fwdModelData);
Az = rocarea(y,truth);

R.y = y;
R.w = w;
R.v = v;
R.y_EEG = y_EEG;
R.fwdModel = fwdModel;
R.dwell = dwell;
R.truth = truth;
R.Az = Az;