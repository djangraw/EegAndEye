% Find fixation duration distribution
%
% Created 3/26/14 by DJ.

% event_types = [R_sq_type(1).regressor_events{end} {'Circle'}];
event_types = {'aD_{0/2}','aT_{0/2}','aD_{1/2}','aT^*_{1/2}','aD_{+/2}','aT_{+/2}','Circle'};
sqnum_types = {'SqNum1','SqNum2','SqNum3','SqNum4','SqNum5'};%,'Circle'};

xDur = 0:.050:1.5;
dur_distribution = zeros(numel(event_types)-1,length(xDur),numel(R_sq_type));
[dur_mean, dur_std] = deal(nan(numel(event_types)-1,numel(R_sq_type)));
for i=1:numel(R_sq_type)
    fprintf('---Subject %d/%d...\n',i,numel(R_sq_type));
    EEG = R_sq_type(i).EEG;
    events = {EEG.event.type};
    latencies = [EEG.event.latency]/EEG.srate;
    trials = [EEG.event.epoch];
    for j=1:numel(event_types)-1
        iThis = find(strcmp(event_types{j},events));        
%         iNot = find(ismember(events,event_types([1:j-1 j+1:end])));
        iSquare = find(ismember(events,sqnum_types));
        duration = nan(1,length(iThis));
        sqnum = nan(1,length(iThis));
        for k=1:numel(iThis)
            if strcmp('Errant',events{iThis(k)-1})
                continue;
            end
            iSqNum = find(latencies(iSquare)==latencies(iThis(k)));
            if isempty(iSqNum)
                continue; 
            end;

            sqnum(k) = find(strcmp(events{iSquare(iSqNum)},sqnum_types));
            if sqnum~=3 %sqnum(k)==1 || sqnum(k)==5
                continue;
            end
%             if isempty(sqnum), break; end;
            iNext_thistrial = iSquare(find(trials(iSquare)==trials(iThis(k)) & iSquare>iThis(k),1));
            if ~isempty(iNext_thistrial)
                duration(k) = latencies(iNext_thistrial(1)) - latencies(iThis(k));
            end
        end
        dur_distribution(j,:,i) = hist(duration,xDur)/sum(~isnan(duration));        
        dur_mean(j,i) = nanmean(duration);
        dur_std(j,i) = nanstd(duration);
        nEvents(j,i) = sum(~isnan(duration));
    end
end

%%
event_types = {'SqNum1','SqNum2','SqNum3','SqNum4','SqNum5','Circle'};
xDur = 0:.010:1.5;
dur_distribution = zeros(numel(event_types)-1,length(xDur),numel(R_sq_type));
for i=1:numel(R_sq_type)
    EEG = R_sq_type(i).EEG;
    events = {EEG.event.type};
    latencies = [EEG.event.latency]/EEG.srate;
    trials = [EEG.event.epoch];
    for j=1:numel(event_types)-1
        iThis = find(strcmp(event_types{j},events));        
        iNext = find(strcmp(event_types{j+1},events));
        duration = nan(1,length(iThis));
        for k=1:numel(iThis)
            iNext_thistrial = iNext(trials(iNext)==trials(iThis(k)));
            if ~isempty(iNext_thistrial)
                duration(k) = latencies(iNext_thistrial(1)) - latencies(iThis(k));
            end
        end
        dur_distribution(j,:,i) = hist(duration,xDur)/sum(~isnan(duration));        
    end
end


%%
figure;
clf; hold on;
% plot(xDur,cumsum(nanmean(dur_distribution,3)'));
Cmap = GetSquaresEventColormap(event_types);
for i=1:size(dur_distribution,1)
    plot(xDur,(nanmean(dur_distribution(i,:,:),3)'),'color',Cmap(i,:),'linewidth',2);
end

% superimpose mean
mean_distrib = mean(cumsum(nanmean(dur_distribution,3)'),2);
hold on;
plot(xDur,mean_distrib,'k','linewidth',2);
% Annotate plot
xlabel('time(seconds)')
ylabel('fraction of trials with duration <= t')
title('Cumulative histogram of event durations')
legend([event_types(1:end-1),{'Mean'}])
ylim([0 1]);
% Estimate median, assuming there's equal numbers of each event type
median_dur = xDur(find(mean_distrib>=0.5,1));
fprintf('median duration is %.3f seconds\n',median_dur);

%%
Cmap = GetSquaresEventColormap(event_types);
Pval_start = nan(numel(event_types)-1);
diffs = nan(numel(event_types)-1);
for i=1:numel(event_types)-1    
    for j=1:numel(event_types)-1
        diffs(i,j) = mean(dur_mean(j,:) - dur_mean(i,:));
        Pval_start(i,j) = signrank(dur_mean(j,:),dur_mean(i,:));
    end
end
% Adjust for FDR
% isHigh = Pval_start>0.5;
% Pval_start(isHigh) = 1-Pval_start(isHigh);        
% Pval = reshape(mafdr(Pval_start(:),'bhfdr',true),size(Pval_start));
% Pval(isHigh) = 1-Pval(isHigh);

p = reshape(mafdr(Pval_start(:),'bhfdr',true),size(Pval_start));

clf;
subplot(1,2,1);
nSubj = numel(R_sq_type);
errorbar(nanmean(dur_mean'),nanstd(dur_mean')/sqrt(nSubj));
set(gca,'xtick',1:numel(event_types)-1,'xticklabel',event_types(1:numel(event_types)-1))
xlabel('Event Type')
ylabel('Fixation Duration (s)');

subplot(1,2,2);
imagesc(diffs);
hold on;
[iStar,jStar] = find(p<0.05);
for i=1:numel(iStar)
    plot(iStar(i),jStar(i),'m*')
end
[iStar,jStar] = find(p<0.005);
for i=1:numel(iStar)
    plot(iStar(i),jStar(i),'y*')
end
xlabel('Event X');
ylabel('Event Y');
title('Dur(Event X) - Dur(Event Y)');
set(gca,'xtick',1:7,'xticklabel',event_types(1:7))
set(gca,'ytick',1:7,'yticklabel',event_types(1:7))
colorbar;

%%
clf;
hold on;
for j=1:numel(event_types)-1
    bar(j,mean(dur_mean(j,:))*1000,'facecolor',Cmap(j,:));
end
errorbar(mean(dur_mean')*1000,std(dur_mean'*1000)/sqrt(nSubj),'k.');
set(gca,'xtick',1:numel(event_types)-1,'xticklabel',event_types(1:numel(event_types)-1))
xlabel('Event Type')
ylabel('Fixation Duration (ms)');
ylim([400 520])

set(gcf,'Position',[1922 535 475 370]);


