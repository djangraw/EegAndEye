function [iObjects_tag_pt, iObjects_eeg_pt,R] = RunResultsThroughTAG(subject,sessions,ALLEEG,cvmode,useica,offset,nSensitivity)

% [iObjects_tag_pt, iObjects_eeg_pt, R] = RunResultsThroughTAG(subject,sessions,ALLEEG,cvmode,useica)
%
% A half-baked script to take the results of the EM jittered logistic
% regression algorithm and see what TAG does with them.  Needs a lot of
% work and cleanup.
%
% INPUTS:
% -subject is the subject number and sessions is a vector of sessions
% numbers of 3DSearch data files as imported with Import_3DS_Data_v3.
% -ALLEEG is a 2-element vector of eeglab structs
% -cvmode is a string
% -useica is a boolean
% -offset is a scalar
% -nSensitivity is a scalar indicating how many times you want to perform
% "self-tuning" (throwing out one predicted target that doesn't match and
% adding another that does). [Default: 0]. fracSensitivity may also be used
% (see RerankObjectsWithTag for details).
%
% Created 5/13/11 by DJ.
% Updated 11/8/11 by DJ - using nTargets/4 eeg predicted targets instead of
% nTargets.
% Updated 8/20/12 by DJ - made a function, made plotting optional
% Updated 3/29/13 by DJ - added cvmode, useica inputs, no 2*stddev cutoff

if nargin<4 || isempty(cvmode)
    cvmode = '10fold';
end
if nargin<5 || isempty(useica)
    useica = false;
end
if nargin<6 || isempty(offset)
    offset = 0; % timing offset (in ms) between EEG and eye tracker
end
if nargin<7 || isempty(nSensitivity)
    nSensitivity = 0;
end

%% Set up
plotresults = 1;

% tWindowStart = 0:100:900;
% tWindowStart = 100:100:900;
tWindowStart = 100:100:1100;
% tWindowStart = 100:50:950;
timesWithOffset = ALLEEG(1).times - offset;
samplesWindowStart = round(interp1(timesWithOffset, 1:ALLEEG(1).pnts, tWindowStart));

tWindowLength = 100;
% tWindowLength = 50;
samplesWindowLength = round(tWindowLength/1000*ALLEEG(1).srate);

[y,w,v,fwdModel] = run_rsvp_classifier(ALLEEG,samplesWindowLength,samplesWindowStart,cvmode,useica);
truth = [zeros(1,ALLEEG(1).trials), ones(1,ALLEEG(2).trials)];
p = y';
Az = rocarea(y,truth);
%% Get results

% load results and params, then...
nTargets = sum(truth);
nDistractors = numel(truth)-nTargets;
trialNumbers = [1:nDistractors, 1:nTargets]; % the trial number (in EEG struct) for each element of p and truth

% sort trials
% [~, best_offset] = max(Azloo); % cherry-pick the window with the best leave-one-out classification
% [p_ordered, order] = sort(p(:,best_offset),1,'descend'); % sort trials in order of descending probability of being a target
[~, order] = sort(p,1,'descend'); % sort trials in order of descending probability of being a target
truth_ordered = truth(order); % sort the truth values the same way

% get eeg-predicted targets by cropping these newly-sorted values
% n_eeg_pt = 8;%round(nTargets/2); % cutoff for being a predicted target (nTargets --> same as target prevalence)
% n_eeg_pt = sum(y>mean(y)+2*std(y))
% if n_eeg_pt<5
%     disp('Using 1.64 std threshold...')
%     n_eeg_pt = sum(y>mean(y)+1.64*std(y))
% end
% if n_eeg_pt<5
%     warning('fewer than 5 trials passed 1.64std threshold!')
    fprintf('Using 1 std threshold... ')
    n_eeg_pt = sum(y>mean(y)+1*std(y));
    fprintf('%d objects passed threshold.\n',n_eeg_pt);
% end
if n_eeg_pt<5
    error('fewer than 5 trials passed 1 std threshold!')    
end
order_eeg_pt = order(1:n_eeg_pt); % crop trial numbers
truth_eeg_pt = truth_ordered(1:n_eeg_pt); % crop truth values

% get eeg struct trial numbers of eeg predicted targets ("eeg_pt")
iTargTrials_eeg_pt = trialNumbers(order_eeg_pt(truth_eeg_pt==1)); % get trial numbers in target EEG struct
iDistTrials_eeg_pt = trialNumbers(order_eeg_pt(truth_eeg_pt~=1)); % get trial numbers in distractor EEG struct

% get a list of the object numbers for the eeg trials
% subject = 2; sessions = 41:48;
[objects, objnames, objlocs, objtimes, objisessions] = GetObjectList(subject,sessions); % get a list of every object seen and its properties
iObjects_targets = EpochToObjectNumber(ALLEEG(2),objtimes, objisessions); % Find the objects that were seen in the target EEG struct's trials
iObjects_distractors = EpochToObjectNumber(ALLEEG(1),objtimes, objisessions); % Find the objects that were seen in the distractor EEG struct's trials

% rerank with TAG
iObjects_eeg_pt = [iObjects_targets(iTargTrials_eeg_pt) iObjects_distractors(iDistTrials_eeg_pt)]; % get object numbers of EEG predicted targets
newRanking = RerankObjectsWithTag(objnames,iObjects_eeg_pt,nSensitivity); % feed these objects as inputs to TAG

% get indices of TAG-predicted targets
n_tag_pt = numel(objects)/4;
iObjects_tag_pt = newRanking(1:n_tag_pt);

%% Display
if plotresults
    % output some stats
    objistarget = strcmp('TargetObject',{objects(:).tag});
    pctCorrect_eeg_pt = sum(objistarget(iObjects_eeg_pt))/numel(iObjects_eeg_pt)*100;
    pctCorrect_tag_pt = sum(objistarget(iObjects_tag_pt))/numel(iObjects_tag_pt)*100;
    pctFound_eeg_pt = sum(objistarget(iObjects_eeg_pt))/sum(objistarget)*100;
    pctFound_tag_pt = sum(objistarget(iObjects_tag_pt))/sum(objistarget)*100;
    fprintf('Percent Correct -- EEG: %.1f, TAG: %.1f\n',pctCorrect_eeg_pt,pctCorrect_tag_pt);
    fprintf('Percent of Targets Found -- EEG: %.1f, TAG: %.1f\n',pctFound_eeg_pt,pctFound_tag_pt);
%     clear pct*
    % Make results struct
    R.pctCorrect_eeg_pt = pctCorrect_eeg_pt;
    R.pctCorrect_tag_pt = pctCorrect_tag_pt;
    R.pctFound_eeg_pt = pctFound_eeg_pt;
    R.pctFound_tag_pt = pctFound_tag_pt;
    R.w = w;
    R.v = v;
    R.y = y;
    R.truth = truth;
    R.fwdModel = fwdModel;
    R.Az = Az;

    figure;
    plot(cumsum(truth_ordered)./(1:numel(truth_ordered)))
    hold on
    PlotVerticalLines([sum(y>mean(y)+1*std(y)), sum(y>mean(y)+1.64*std(y)), sum(y>mean(y)+2*std(y))],'m:');
    plot(get(gca,'xlim'),[0.25 0.25],'k--');
    xlabel('# trials included');
    ylabel('Precision of set');
    title(sprintf('S%d EEG precision',subject));

%     %% Run ImageAllSessions on results
%     levelname = 'GridHuge.png'; % 'SnakeHuge.png'
%     levelinfo=load(levelname(1:end-4));
%     figure(356); clf
%     ImageAllSessions(subject,sessions,levelname,levelinfo.zeroPoint,iObjects_eeg_pt);
end