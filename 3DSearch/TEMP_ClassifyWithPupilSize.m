% Created 7/9/13 by DJ.
% Updated 7/11/13 by DJ - fixed bug from subjects w/ fewer sessions
%% Set up
subjects = [22:30 32];
sessions_cell = {2:14, [3 6:17], 1:15, 1:15, 1:15, 1:15, 1:15, 1:15, [1:10 12:15], 2:16};
offsets = [-12 -4 48 60 68 68 92 112 0 0];
epochRange = [-1000 4000];

xeye_all = cell(1,numel(subjects));
Az_ps_all = [];
ps_median_targ_all = [];
ps_median_dist_all = [];

for iSubj = 1:numel(subjects);
    subject = subjects(iSubj);
    sessions = sessions_cell{iSubj};
    offset = offsets(iSubj);
    
    %% Load EEG and behavior structs
    y = loadBehaviorData(subject,sessions,'3DS');
    EEG = pop_loadset('filename',...
        sprintf('3DS-%d-all-filtered-noduds-noeog-epoched-ica.set',subject));

    %% Get eye pos and pupil size
    [ps,xeye,yeye] = deal(cell(1,numel(sessions)));
    for i=1:numel(sessions)        
        load(sprintf('3DS-%d-%d-eyepos',subject,sessions(i)));
        xeye{i} = eyepos(:,1);
        yeye{i} = eyepos(:,2);
        ps{i} = InterpolateBlinks(pupilsize,y(i).eyelink.record_time-1+(1:length(pupilsize)),y(i));
    end 
    
    %% Get epochs
    [ps_epoch,t_epoch,isTargetEpoch] = GetEpochedPupilSize(y,ps,EEG,epochRange);
    xeye_epoch = GetEpochedPupilSize(y,xeye,EEG,epochRange);

    %% Normalize pupil size to be percentage change
    mean_ps = nanmean(cat(1,ps{:}));
    ps_epoch_pct = ps_epoch/mean_ps*100;

    %% Remove baseline
    tBaseline = [-1000 0];%[0 100];
    isBase = t_epoch>=tBaseline(1) & t_epoch<=tBaseline(2);
    baseline = nanmean(ps_epoch_pct(:,isBase),2);
    ps_epoch_pct_baserem = ps_epoch_pct - repmat(baseline,1,size(ps_epoch,2));

    %% Get AUC values
    Az_ps = nan(1,size(ps_epoch,2));
    for i=1:size(ps_epoch,2)
        Az_ps(i) = rocarea(ps_epoch_pct_baserem(:,i),isTargetEpoch);
    end

    %% Plot ps
    clims = [min(ps_epoch_pct_baserem(:)), max(ps_epoch_pct_baserem(:))];
    figure;

    subplot(4,1,1); hold on
    imagesc(t_epoch,1:sum(isTargetEpoch),ps_epoch_pct_baserem(isTargetEpoch==1,:));
    plot([0 0],get(gca,'ylim'),'k-');
    set(gca,'clim',clims);
    axis([t_epoch(1), t_epoch(end),0.5,sum(isTargetEpoch)+0.5]);
    colorbar
    ylabel(sprintf('pupil area change\n(%% of subj mean)'))
    title('Target Epochs')

    subplot(4,1,2); hold on
    imagesc(t_epoch,1:sum(~isTargetEpoch),ps_epoch_pct_baserem(isTargetEpoch==0,:));
    plot([0 0],get(gca,'ylim'),'k-');
    set(gca,'clim',clims)
    axis([t_epoch(1), t_epoch(end),0.5,sum(~isTargetEpoch)+0.5]);
    colorbar
    ylabel(sprintf('pupil area change\n(%% of subj mean)'))
    title('Distractor Epochs')

    subplot(4,1,3); hold on;
    ps_median_targ = nanmedian(ps_epoch_pct_baserem(isTargetEpoch==1,:),1);
    ps_stderr_targ = nanstd(ps_epoch_pct_baserem(isTargetEpoch==1,:),[],1)/sqrt(sum(isTargetEpoch==1));
    ps_median_dist = nanmedian(ps_epoch_pct_baserem(isTargetEpoch==0,:),1);
    ps_stderr_dist = nanstd(ps_epoch_pct_baserem(isTargetEpoch==0,:),[],1)/sqrt(sum(isTargetEpoch==0));

    % plot
    % plot(t_epoch,[ps_median_targ; ps_median_dist; ps_median_targ - ps_median_dist]);
    plot(t_epoch,ps_median_targ,'r','linewidth',2);
    plot(t_epoch,ps_median_dist,'b','linewidth',2);
    plot(t_epoch,ps_median_targ - ps_median_dist,'g','linewidth',2);
    ErrorPatch(t_epoch,ps_median_targ,ps_stderr_targ,'r','r');
    ErrorPatch(t_epoch,ps_median_dist,ps_stderr_dist,'b','b');
    plot(t_epoch,ps_median_targ - ps_median_dist,'g','linewidth',2);
    plot([t_epoch(1) t_epoch(end)],[0 0],'k-');
    plot([0 0],get(gca,'ylim'),'k-');
    set(colorbar,'Visible','off') % dummy colorbar to make axes line up
    legend('Target (+/- stderr)','Distractor(+/- stderr)','Difference','Location','NorthWest');
    ylabel(sprintf('Median pupil area change\n(%% of subj mean)'))

    subplot(4,1,4); hold on
    plot(t_epoch,Az_ps);
    set(colorbar,'Visible','off') % dummy colorbar to make axes line up
    ylabel('AUC')
    xlabel('Time from saccade to object (ms)');
    plot([t_epoch(1) t_epoch(end)],[0.5 0.5],'k-');
    plot([0 0],get(gca,'ylim'),'k-');

    MakeFigureTitle(sprintf('Subject %d Pupil Area ([%d %d] ms baseline)',subject,tBaseline(1),tBaseline(2)),0);
    
    % Add to results
    xeye_all{iSubj} = xeye_epoch;
    Az_ps_all(iSubj,:) = Az_ps;
    ps_median_targ_all(iSubj,:) = ps_median_targ;
    ps_median_dist_all(iSubj,:) = ps_median_dist;
end

%% Plot cross-subject averages
figure;
subplot(2,1,1); hold on;
mtarg = mean(ps_median_targ_all,1);
mdist = mean(ps_median_dist_all,1);
mdiff = mean(ps_median_targ_all - ps_median_dist_all,1);
starg = std(ps_median_targ_all,[],1)/sqrt(numel(subjects));
sdist = std(ps_median_dist_all,[],1)/sqrt(numel(subjects));
sdiff = std(ps_median_targ_all - ps_median_dist_all,[],1)/sqrt(numel(subjects));


plot(t_epoch,mtarg,'r','linewidth',2);
plot(t_epoch,mdist,'b','linewidth',2);
plot(t_epoch,mdiff,'g','linewidth',2);
ErrorPatch(t_epoch,mtarg,starg,'r','r');
ErrorPatch(t_epoch,mdist,sdist,'b','b');
ErrorPatch(t_epoch,mdiff,sdiff,'g','g');
plot([t_epoch(1) t_epoch(end)],[0 0],'k-');
plot([0 0],get(gca,'ylim'),'k-');
set(colorbar,'Visible','off') % dummy colorbar to make axes line up
legend('Target','Distractor','Difference','Location','NorthWest');
ylabel(sprintf('Median pupil area change\n(%% of subj mean)\nmean +/- stderr across subjects'))

subplot(2,1,2); hold on
ErrorPatch(t_epoch,mean(Az_ps_all,1),std(Az_ps_all,[],1)/sqrt(numel(subjects)),'b','b');
set(colorbar,'Visible','off') % dummy colorbar to make axes line up
ylabel(sprintf('AUC\nmean +/- stderr across subjects'))
xlabel('Time from saccade to object (ms)');
plot([t_epoch(1) t_epoch(end)],[0.5 0.5],'k-');
plot([0 0],get(gca,'ylim'),'k-');
