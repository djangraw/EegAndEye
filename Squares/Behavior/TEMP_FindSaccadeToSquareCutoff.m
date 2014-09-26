% TEMP_FindSaccadeToSquareCutoff.m
%
% Find the cutoff that best captures the subject's saccade choices when
% makeing multiple saccades to a square: that is, that balances the number
% with >1 saccade above the cutoff and with 0 saccades above the cutoff.
% 
% Created 2/13/13 by DJ for one-time use.

% Set up
subjects = [9:11 13:15 17:19];
clear data
for i=1:numel(subjects)
    foo = load(sprintf('sq-%d-all',subjects(i)));
    data(i,:) = foo.y;
end
%% Run FindSaccadesToSquares
cutoffs = 0:10:500;
clear nones multis
for i=1:size(data,1) % subject
    for j=1:size(data,2) % session
        for k=1:numel(cutoffs) % cutoff
            [~,~,nones(i,j,k),multis(i,j,k)] = FindSaccadesToSquares(data(i,j),cutoffs(k));
        end
    end
end

%% Plot mean results
% Get means
mean_nones = squeeze(mean(nones,2));
mean_multis = squeeze(mean(multis,2));
meanmean_nones = mean(mean_nones,1);
meanmean_multis = mean(mean_multis,1);
% Plot
cla; hold on
plot(cutoffs, mean_nones,'-');
plot(cutoffs, mean_multis,'--');
plot(cutoffs, meanmean_nones, 'k-', 'linewidth',2);
plot(cutoffs, meanmean_multis, 'k--', 'linewidth',2);

legendText = cell(1,numel(subjects));
for i=1:numel(subjects)
    legendText{i} = sprintf('S%d',subjects(i));
end
legend(legendText)
xlabel('fixation duration cutoff (ms)');
ylabel('# trials across all sessions');
title(sprintf('Saccade-To-Square fixation duration cutoff\n(solid: squares w/ none exceeding cutoff. dotted: squares w/ multiple exceeding cutoff)\n(Black: mean across subjects)'));


%% find intersection point
iOver = find(meanmean_nones>meanmean_multis,1);
iUnder = iOver-1;
% Solve linear equation
m_nones = (meanmean_nones(iOver)-meanmean_nones(iUnder))/(cutoffs(iOver)-cutoffs(iUnder));
m_multis = (meanmean_multis(iOver)-meanmean_multis(iUnder))/(cutoffs(iOver)-cutoffs(iUnder));
b_nones = meanmean_nones(iOver) - m_nones*cutoffs(iOver);
b_multis = meanmean_multis(iOver) - m_multis*cutoffs(iOver);
% Calculate intersection
cutoffs_intersect = (b_nones-b_multis)/(m_multis-m_nones);

% plot green star at intersection
plot(cutoffs_intersect,m_multis*cutoffs_intersect+b_multis,'g*');
% display
fprintf('Cutoff intersection: %.2f\n',cutoffs_intersect)



