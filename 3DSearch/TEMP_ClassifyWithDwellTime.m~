% TEMP_ClassifyWithDwellTime
%
% Created 5/22/13 by DJ for one-time use.

%% Set up

subjects = [22:30 32];
sessions_cell = {2:14, [3 6:17], 1:15, 1:15, 1:15, 1:15, 1:15, 1:15, [1:10 12:15], 2:16};
offsets = [-12 -4 48 60 68 68 92 112 0 0];
cvmode = '10fold';

%% Classify with dwell times alone
Az_dwell = zeros(1,numel(subjects));
% dwell = cell(1,numel(subjects));
% isToTarget = cell(1,numel(subjects));
[y_dwell,truth_dwell] = deal(cell(1,numel(subjects)));
for i=1:numel(subjects)
    fprintf('--- Subject %d/%d: ---\n',i,numel(subjects));
    [dwell{i} isToTarget{i}] = GetDwellTimes(subjects(i),sessions_cell{i});
    dwellz = [dwell{i}{:}];
    isTargetz = [isToTarget{i}{:}];
    y_dwell{i} = dwellz(~isnan(dwellz));
    truth_dwell{i} = isTargetz(~isnan(dwellz));
%     Az_dwell(i) = rocarea(y_dwell{i},truth_dwell{i});
end
    

%% Try classifying with both
clear R;
iObjects_eeg_pt = cell(1,numel(subjects));
for i=10%1:numel(subjects)
    fprintf('--- Subject %d/%d ---\n',i,numel(subjects));
%     load(sprintf('ALLEEG%d_eyeposcorrected.mat',subjects(i))); % ALLEEG
%     load(sprintf('ALLEEG%d_EpochedIcaCropped.mat',subjects(i))); % ALLEEG
    load(sprintf('ALLEEG%d_NoeogEpochedIcaCropped.mat',subjects(i))); % ALLEEG
    y = loadBehaviorData(subjects(i),sessions_cell{i},'3DS');
    % --- IF DWELL TIME IS KNOWN --- %
%     R(i) = ClassifyWithEegAndDwellTime(ALLEEG,y,cvmode,offsets(i),...
%       R_both_eica_nested(i).dwell); % if dwell time is known
    % --- IF DWELL TIME IS UNKNOWN --- %
%     R(i) = ClassifyWithEegAndDwellTime(ALLEEG,y,cvmode,offsets(i)); % if dwell time is unknown
    % --- IF DWELL TIME IS NOT USED --- %
    R(i) = ClassifyWithEegAndDwellTime(ALLEEG,y,cvmode,offsets(i),[]); % for no dwell times
    iObjects_eeg_pt{i} = GetEegPredictedTargets(subjects(i),sessions_cell{i},ALLEEG,R(i).y,R(i).truth);
end
disp('---DONE!---');

%% Compare results of dwell time vs. eeg classification
Az_dwell = Az_dwell;
Az_EEG = [R_9bin.Az];
Az_both = [R_both.Az];
Az_eeg_top10 = [R_9bin_top10.Az];
Az_both_top10 = [R_both_top10.Az];
Az_eeg_nested = [R_EEG_nested.Az];
Az_both_nested = [R_both_nested.Az];
Az_eeg_top10_nested = [R_EEG_top10_nested.Az];
Az_both_top10_nested = [R_both_top10_nested.Az];
Az_both_nested_20fold = [R_both_nested_20fold.Az];
Az_both_eica_nested = [R_both_eica_nested.Az];
Az_both_eica_top10_nested = [R_both_eica_top10_nested.Az];
Az_both_noeog_eica20_nested = [R_both_noeog_eica20_nested.Az];
Az_eeg_noeog_eica20_nested = [R_EEG_noeog_eica20_nested.Az];
% Az_new = [R.Az];

% Az_all  = [Az_dwell; Az_EEG; Az_both; Az_eeg_top10; Az_both_top10; Az_eeg_nested; Az_both_nested; Az_eeg_top10_nested; Az_both_top10_nested; Az_both_nested_20fold; Az_both_eica_nested; Az_both_eica_top10_nested; Az_eeg_noeog_eica20_nested; Az_both_noeog_eica20_nested];
% Az_all  = [Az_dwell; Az_eeg_nested; Az_both_nested; Az_both_nested_20fold];
Az_all = [Az_dwell; Az_eeg_noeog_eica20_nested; Az_both_noeog_eica20_nested];

cla; hold on
% plot([Az_all]','.-')
% Plot normally
% PlotUniqueLines(1:size(Az_all,2),[Az_all]','.')
% Plot ordered by EEG Az
[~,order] = sort(Az_eeg_noeog_eica20_nested,'descend');
PlotUniqueLines(1:size(Az_all,2),Az_all(:,order)','.');
% Annotate plot
plot([1 numel(subjects)],[0.5 0.5],'k--')
% xlabel('subject')
xlabel('subject (sorted by EEG Az score)')
ylabel('Az')
title('Eye position vs. EEG')
ylim([0.3 1])

% legend('dwell time','EEG (9-bin HDCA)','dwell + EEG','EEG (top 10 ICs)',...
%   'dwell + EEG (top 10 ICs)','EEG (nested 10-fold)','both (nested 10-fold)','EEG (top 10, nested)','both (top 10, nested)','both (nested 20-fold)','both (epoched ICA, nested)','both (eICA, top 10, nested)','EEG (no eog eICA20, nested)','both (no eog eICA20, nested)');
% legend('dwell time','EEG (nested 10-fold)','both (nested 10-fold)','both (nested 20-fold)');
legend('dwell time','EEG (no eog epoched ICA-20, nested)','both (no eog eICA20, nested)');

%% Fancy:

% Get data
load('TEMP_EegAndDwell_Results.mat');
Az_both_noeog_eica20_nested = [R_both_noeog_eica20_nested2.Az];
Az_eeg_noeog_eica20_nested = [R_EEG_noeog_eica20_nested2.Az];
Az_all = [Az_dwell; Az_eeg_noeog_eica20_nested; Az_both_noeog_eica20_nested];

% Set plot parameters
fontname = 'Futura'; % Futura, Helvetica, Gotham are all good
fontsize = 15;
linewidth = 2;
markersize = 20;
nSubjects = numel(subjects);

% Set up
cla; hold on;
set(gca,'xtick',1:nSubjects,'ytick',0.3:.1:1,'box','on','fontname',fontname,'fontsize',fontsize)
% Plot
PlotUniqueLines(1:nSubjects, Az_all', '.', linewidth, markersize)
plot([0 nSubjects+1],[0.5 0.5],'k--','LineWidth',linewidth)
% Annotate
xlabel('subject')
ylabel('Area Under ROC Curve')
title('Classifier Performance')
ylim([0.3 1])
legend('Dwell Time','EEG','Hybrid','Chance')