% TEMP_QuantifyTimingOffsets.m
%
% Created 5/15/13 by DJ.
% Updated 7/23/13 by DJ - added loading code, 2 subjects.

%% Load datasets
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
subjects = [22:30 32];
nSubjects = numel(subjects);
for i=1:nSubjects
    EEG = pop_loadset('filename',sprintf('3DS-%d-all-filtered-noduds.set',subjects(i)));
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
end

%% Get ERPs for saccade-to-distractor trials (include eye components!)
[OzErp,chanlabels,ERP] = deal(cell(1,nSubjects));
for i=1:nSubjects
    fprintf('S%d/%d...\n',i,nSubjects)
    EEG2 = pop_epoch(ALLEEG(i),{'410'},[-.5 1]);
    chanlabels{i} = {EEG2.chanlocs.labels};
    ERP{i} = mean(EEG2.data,3);
    iOz = strcmp('O1',chanlabels{i});
    OzErp{i} = mean(EEG2.data(iOz,:,:),3);
end
times = EEG2.times;
%% Get ERPs on matching electrodes
chans_all = chanlabels{1};
for i=1:nSubjects
    chans_all = intersect(chanlabels{i},chans_all);
end
ERPmatch = nan(length(chans_all),size(ERP{1},2),nSubjects);
for i=1:nSubjects
    ERPmatch(:,:,i) = ERP{i}(ismember(chanlabels{i},chans_all),:);
end
chanlocs = ALLEEG(1).chanlocs(ismember(chanlabels{1},chans_all));

%% Plot matched ERPs as scalp maps
legendstr = cell(1,nSubjects);
for i=1:nSubjects
    legendstr{i} = sprintf('S%d',subjects(i));
end

figure(45);
bintimes = -50:100:500;
results = GetScalpMaps(ERPmatch,times,bintimes,100);
PlotScalpMaps(results,chanlocs,[],bintimes,legendstr);

%% Plot topo movie
for i=1:nSubjects
    MakeTopoMovie(ERP{i},times,ALLEEG(i).chanlocs,[]);
    MakeFigureTitle(sprintf('S%d',subjects(i)),0);
end
%% Plot ERPs on representative electrode 
figure(46); clf;
plot(times,cat(1,OzErp{:})')
xlim([-100 500])
legend(legendstr);
%% Plot ERP mean over all electrodes
figure(47); cla; hold on;
colors = {'b','r','g','c','m','y','k','b:'};    
for i=1:nSubjects
    meanERP(i,:) = mean(ERP{i},1);
    plot(times,meanERP(i,:),colors{i});
end
xlim([-100 500])
legend(legendstr);
%% Use Woody algorithm to line up trials
baselineWin = [100 0];
finalMatchStrength = LineUpTrials_rawdata(ERPmatch,times,[0 500],10,chanlocs);

%% Plot results
[~,iMax] = max(finalMatchStrength,[],2);

for i=1:nSubjects
    plot(EEG2.times(iMax(i)),meanERP(i,iMax(i)),[colors{i} '.']);
end