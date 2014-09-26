function [all_dist,all_targ,times,chanlocs_all] = TEMP_Make3dsErps(subjects,sessions_cell,offsets,gridOfChans)

% Created 7/16/13 by DJ.
% Updated 7/18/13 by DJ - fixed offset alignment
% Updated 7/23/13 by DJ - added outputs

nSubj = numel(subjects);
for i=1:nSubj
    fprintf('Loading data from S%d/%d...\n',i,nSubj);
%     y = loadBehaviorData(subjects(i),sessions_cell{i},'3DS');
    load(sprintf('ALLEEG%d_NoeogEpochedIcaCropped',subjects(i)));
    dist{i} = mean(ALLEEG(1).data,3);
    targ{i} = mean(ALLEEG(2).data,3);
    diffs{i} = targ{i}-dist{i};
    t_erp{i} = ALLEEG(1).times-offsets(i);
    chanlocs{i} = ALLEEG(1).chanlocs;
end

disp('Calculating...')
% Reconcile chanlocs
channies = [chanlocs{:}];
chans = unique({channies.labels});
for i=1:numel(subjects)
    chans = chans(ismember(chans,{chanlocs{i}.labels}));
end
chanlocs_all = chanlocs{1}(ismember({chanlocs{1}.labels},chans));


for i=1:nSubj
    dist_cropped{i} = dist{i}(ismember({chanlocs{i}.labels},chans),:);
    targ_cropped{i} = targ{i}(ismember({chanlocs{i}.labels},chans),:);
    diffs_cropped{i} = diffs{i}(ismember({chanlocs{i}.labels},chans),:);
end

all_dist = cat(4,dist_cropped{:});
all_targ = cat(4,targ_cropped{:});
all_diffs = cat(4,diffs_cropped{:});
% mean_dist = mean(cat(3,dist_cropped{:}),3);
% std_dist = std(cat(3,dist_cropped{:}),[],3);
% mean_targ = mean(cat(3,targ_cropped{:}),3);
% std_targ = std(cat(3,targ_cropped{:}),[],3);
% mean_diff = mean(cat(3,diffs_cropped{:}),3);
% std_diff = std(cat(3,diffs_cropped{:}),[],3);

disp('Plotting...')
% Set up figure
% figure;
clf;
% Plot
% PlotResponseFnsGrid(cat(3,all_dist,all_targ,all_diffs),{'Dist', 'Targ','Targ - Dist'},t_erp{1},chanlocs_all,gridOfChans);
times = t_erp{1};
PlotResponseFnsGrid(cat(3,all_dist,all_targ),{'Distractors', 'Targets'},times,chanlocs_all,gridOfChans);

%% Annotate plot
% Set plot parameters
fontname = 'Futura';
fontsize = 15;
% linewidth = 2;
% markersize = 20;
gridOfChansT = gridOfChans';
for i=1:numel(gridOfChans)
    subplot(size(gridOfChans,1),size(gridOfChans,2),i);
    set(gca,'box','on','fontname',fontname,'fontsize',fontsize);
    xlabel('Time from First Saccade to Object (ms)')
    ylabel('Response (\muV)')
    title(gridOfChansT{i});
end
disp('Done!');