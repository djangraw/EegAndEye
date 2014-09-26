% RunPupilResultsThroughTAG.m
%
% A half-baked script to take the results of the EM jittered logistic
% regression algorithm and see what TAG does with them.  Needs a lot of
% work and cleanup.
%
% Created 11/7/11 by DJ.


%% Get results
% subject = 18; sessions = 1:10;
% Get pupil responses
h = PlotEyeErps_MultiSession(subject,sessions,1,[-500 2000]);
tp = h.data.targetPupEpochs;
dp = h.data.distractorPupEpochs;
% Subtract baseline
isBaseline = (h.data.epochTimes>-200 & h.data.epochTimes<=0);
for i=1:size(tp,1)
    tp(i,:) = tp(i,:) - mean(tp(i,isBaseline));
end
for i=1:size(dp,1)
    dp(i,:) = dp(i,:) - mean(dp(i,isBaseline));
end
% find best window for discrimination
for i=1:size(tp,2)
    Az(i) = rocarea([tp(:,i); dp(:,i)], [ones(size(tp,1),1); zeros(size(dp,1),1)]);
end

%% Extract important info
% Get y and truth info
y = [dp; tp];
truth = [zeros(size(dp,1),1); ones(size(tp,1),1)];
eventSessions = [h.data.distractorEventSessions; h.data.targetEventSessions];
eventTimes = [h.data.distractorEventTimes; h.data.targetEventTimes];

% sort trials
[~, best_offset] = max(Az); % cherry-pick the window with the best leave-one-out classification
[y_ordered, order] = sort(y(:,best_offset),1,'descend'); % sort trials in order of descending probability of being a target
% Get top results
top = order(1:round(end/4)); % label top 25% as EEG-predicted targets
y_top = y(top,best_offset);
truth_top = truth(top); % sort the truth values the same way
eventSessions_top = eventSessions(top);
eventTimes_top = eventTimes(top);

% Get object info
[objects, objnames, objlocs, objtimes, objsessions] = GetObjectList(subject,sessions,'eyelink'); % get a list of every object seen and its properties
iObjects_eeg_pt = zeros(size(top)); % eeg predicted targets
for i=1:numel(top)
    iObj = find(objsessions==eventSessions_top(i) & objtimes==eventTimes_top(i)/1000);    
    if ~isempty(iObj)
        iObjects_eeg_pt(i) = iObj;
    else
        warning('Object not found!')
    end
end

%% Rerank with TAG
newRanking = RerankObjectsWithTag(objnames,iObjects_eeg_pt); % feed these objects as inputs to TAG

% get indices of TAG-predicted targets
n_tag_pt = round(numel(objects)/4);
iObjects_tag_pt = newRanking(1:n_tag_pt);

%% Display

% output some stats
objistarget = strcmp('TargetObject',{objects(:).tag});
pctCorrect_eeg_pt = sum(objistarget(iObjects_eeg_pt))/numel(iObjects_eeg_pt)*100;
pctCorrect_tag_pt = sum(objistarget(iObjects_tag_pt))/numel(iObjects_tag_pt)*100;
pctFound_eeg_pt = sum(objistarget(iObjects_eeg_pt))/sum(objistarget)*100;
pctFound_tag_pt = sum(objistarget(iObjects_tag_pt))/sum(objistarget)*100;
fprintf('Percent Correct -- EEG: %.1f, TAG: %.1f\n',pctCorrect_eeg_pt,pctCorrect_tag_pt);
fprintf('Percent of Targets Found -- EEG: %.1f, TAG: %.1f\n',pctFound_eeg_pt,pctFound_tag_pt);
clear pct*

%% Run ImageAllSessions on results
clf;
ImageAllSessions(subject,sessions,'GridHuge.png',[15,9.5],iObjects_eeg_pt);