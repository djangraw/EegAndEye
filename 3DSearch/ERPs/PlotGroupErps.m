% PlotGroupErps
% 
% Created 12/16/10 by DJ.
% Updated 1/7/11 by DJ - playing with constants, mostly
% Updated 7/27/11 by DJ - work with new version of RemoveEyeBlinkTrials and
% ExcludeTrialsByCutoff.

%% EXTRACT
subjects = [6 7 2];

for i=1:numel(subjects)
    subject = subjects(i);
    EpochEeglabFile;
    GetNumbers;
    for j=2:5 % the sets that will be saved later
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',j,'study',0); 
        EEG = RemoveEyeBlinkTrials(EEG,Numbers.BLINK,[-350 1000],[]);
%         PlotExclusionHistogram(EEG.data,EEG.times,[0 1000]);
%         EEG = ExcludeTrialsByCutoff(EEG,100,[-350 1000],[]);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    end
    
    % Save results
    pop_saveset(ALLEEG(2), 'filename',sprintf('%d-TargetSac.set',subject));
    pop_saveset(ALLEEG(3), 'filename',sprintf('%d-DistractorSac.set',subject));
    pop_saveset(ALLEEG(4), 'filename',sprintf('%d-TargetApp.set',subject));
    pop_saveset(ALLEEG(5), 'filename',sprintf('%d-DistractorApp.set',subject));
end

%% LOAD 
% Clear datasets from EEGLAB
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
eeglab redraw;

for i=1:numel(subjects)
    subject = subjects(i);
    % Load results    
    EEG = pop_loadset('filename',sprintf('%d-TargetApp.set',subject));
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = pop_loadset('filename',sprintf('%d-DistractorApp.set',subject));
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
end

eeglab redraw;

%% PLOT ALL CHANNEL ERPS

% Plot ERPs
[erp1 erp2 erpsub time] = pop_comperp( ALLEEG, 1,1:2:(2*numel(subjects)),2:2:(2*numel(subjects)),'addavg','on','subavg','on','diffavg','off','diffstd','off','lowpass',30,'tplotopt',{'ydir',1});

% Make TopoMovie figure
MakeTopoMovie(erpsub,time,ALLEEG(1).chanlocs,0.020);


%% PLOT CHANNEL ERPS
channels = [6 15 25]; % Fz Cz Pz, for 'noduds' datasets
% channels = 1:74;
allSigTimes = [];
for i = 1:numel(channels)
    dataTarg = [];
    dataDis = [];
    for j=1:numel(subjects)
        % For a trial-level average...
        dataTarg = [dataTarg; permute(ALLEEG(2*j-1).data(channels(i),:,:),[3 2 1])];
        dataDis = [dataDis; permute(ALLEEG(2*j).data(channels(i),:,:),[3 2 1])];
        
        % For a grand average...
    %     dataTarg(j,:) = mean(ALLEEG(2*j-1).data(channel,:,:),3); 
    %     dataDis(j,:) = mean(ALLEEG(2*j).data(channel,:,:),3);
    end
    % Channel ERPs
    figure; 
    [~, ~, sigTimes] = PlotChannelErp(dataTarg,dataDis,ALLEEG(1).times,ALLEEG(1).chanlocs(channels(i)).labels);
    % Resize figure
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1:2) 600 250]);
%     xlim('auto')
%     close;
%     allSigTimes = [allSigTimes sigTimes];
end
clear channels dataTarg dataDis i j

%% SCALP MAPS

% timepoints = [.050 .150 .250 .350 .450];
% timepoints = [.550 .650 .750 .850 .950];
timepoints = [ .150 .450 .750];
color_limits = [-8 8];
figure;
nRows = numel(timepoints);
set(gcf,'Position',[0 200 550 200*nRows])
for i=1:numel(timepoints)
    subplot(nRows,2,2*i-1);
    iTime = find(time>=timepoints(i),1);
    topoplot(erp1(:,iTime),ALLEEG(1).chanlocs,'MapLimits',color_limits,'conv','on');
    title(sprintf('Targets, t = %g ms',time(iTime)*1000));  
    subplot(nRows,2,2*i);
    topoplot(erp2(:,iTime),ALLEEG(1).chanlocs,'MapLimits',color_limits,'conv','on');
    title(sprintf('Distractors, t = %g ms',time(iTime)*1000));  
end
colorbar;
clear timepoints color_limits nRows iTime i
