% Process3dsData_MultiSubject.m
%
% Created 5/13/13 by DJ
% Updated 5/16/13 by DJ - added offset, fracSensitivity

subjects = [22:30 32];
folders = {'2013-03-27-Pilot-S22', '2013-03-28-Pilot-S23', ...
    '2013-04-29-Pilot-S24', '2013-05-01-Pilot-S25', '2013-05-02-Pilot-S26',...
    '2013-05-03-Pilot-S27', '2013-05-07-Pilot-S28', '2013-05-10-Pilot-S29',...
    '2013-05-15-Pilot-S30', '2013-06-06-Pilot-S32'};
sessions_cell = {2:14, [3 6:17], 1:15, 1:15, 1:15, 1:15, 1:15, 1:15, [1:10 12:15], 2:16};
offsets = [-12 -4 48 60 68 68 92 112 0 0]; % CHECK S30 and S32
fracSensitivity = 0.9999999; % used to calculate nSensitivity
homedir = '/Users/dave/Documents/Data/3DSearch';
cvmode = '10fold';

% set params
pixelThreshold = 100;
timeLimits = [0 Inf];

%% Apply eye pos correction
cd(homedir)
for i=10%1:numel(subjects)   
    cd(folders{i});
    subject = subjects(i);
    % get session #s
    sessions = sessions_cell{i};
    
    % check current calibration
    clear y
    for j=1:numel(sessions)
        foo = load(sprintf('3DS-%d-%d.mat',subject,sessions(j)));
        y(j) = foo.x;
    end
    calz = [y.calibration];
    figure(999); cla; hold on;
    plot([calz.eye_offset_x],[calz.eye_offset_y],'r.');
    drawnow;
    
    % Do correction
    calibration = TestEyePosCorrection(subject,sessions,pixelThreshold);
    figure(999);
    plot([calibration.eye_offset_x],[calibration.eye_offset_y],'.')
    title(sprintf('S%d',subject)); xlabel('x offset'); ylabel('y offset'); legend('current','new');
    drawnow;
    
    foo = input('Apply Correction (Y/N/NULL)?>>','s');
    if strcmp(foo,'Y')
        disp('CHANGING CALIBRATION!') 
        ApplyEyePosCorrection(subject,sessions,calibration,pixelThreshold,timeLimits);    
    elseif strcmp(foo,'NULL')
        disp('APPLYING NULL CALIBRATION!')  
        ApplyEyePosCorrection(subject,sessions,[],pixelThreshold,timeLimits); 
    else
        disp('NOT CHANGING CALIBRATION!')        
    end
    cd ..
end

%% Run Processing
clear R iObjects_eeg_pt
cd(homedir)
for i=1:numel(subjects)
    cd(folders{i});
    [R(i), iObjects_eeg_pt{i}] = Process3dsData(subjects(i),sessions_cell{i},offsets(i),fracSensitivity,cvmode);
    cd ..
end

%% Get TSP and stats
clear stats
usegridconstraints = true;
levelname = 'GridHuge.jpg';
for i=1:numel(subjects)
    fprintf('--- Subject %d/%d ---\n',i,numel(subjects));
    subject = subjects(i);    
    sessions = sessions_cell{i};
    stats(i) = CalculateImprovement(subject,sessions,levelname,iObjects_eeg_pt{i},usegridconstraints,fracSensitivity);
end

%% SET CHANCE LEVELS

chance.precision = [25 25];
chance.pctFound = [20/56*.25*100 25];
chance.pctDistance = 100;

%% PLOT
% Set plot parameters
fontname = 'Futura';
fontsize = 15;
linewidth = 2;
markersize = 20;

subjectstrings = cell(1,numel(subjects));
for i=1:numel(subjects)
    subjectstrings{i} = sprintf('S%d',i);
end

% Make plots
figure(251); clf;
subplot(2,3,1); cla; hold on;
% plot(cat(1,stats.precision)', '.-','linewidth',linewidth,'markersize',markersize);
PlotUniqueLines([1 2], cat(1,stats.precision)','.',linewidth,markersize);
plot(chance.precision,'k.--','linewidth',linewidth,'markersize',markersize);
set(gca,'xtick',[1 2],'xticklabel',{'EEG','EEG-TAG'},'box','on','fontname',fontname,'fontsize',fontsize)
xlim([0 3]);
ylim([0 100]);
ylabel('% Precision');
xlabel('Target Prediction System');
title('Precision of Predicted Targets');
legend([subjectstrings, {'Chance'}],'Location','SouthEast')

% figure(252); clf;
subplot(2,3,2); cla; hold on;
% plot(cat(1,stats.pctFound)', '.-','linewidth',linewidth,'markersize',markersize);
PlotUniqueLines([1,2],cat(1,stats.pctFound)','.',linewidth,markersize);
plot(chance.pctFound,'k.--','linewidth',linewidth,'markersize',markersize);
% set(gca,'xtick',[1 2],'xticklabel',{'Training','TAG'},'box','on','fontname',fontname,'fontsize',fontsize)
set(gca,'xtick',[1 2],'xticklabel',{'EEG','EEG-TAG'},'box','on','fontname',fontname,'fontsize',fontsize)
xlim([0 3]);
ylim([0 100]);
xlabel('Target Prediction System');
ylabel('% Targets Identified')
title('Percentage of True Targets Identified');
% legend([subjectstrings, {'Chance'}],'Location','SouthEast')

% figure(253); clf;
subplot(2,3,3); cla; hold on;
pctFound = cat(1,stats.pctFound);
% plot([zeros(numel(stats),1), cat(1,stats.pctDistance)]', [zeros(numel(stats),1), pctFound(:,2)]', '.-','linewidth',linewidth,'markersize',markersize);
PlotUniqueLines([zeros(numel(stats),1), cat(1,stats.pctDistance)]', [zeros(numel(stats),1), pctFound(:,2)]','.',linewidth,markersize);
plot([0 100], [0 chance.pctDistance],'k.--','linewidth',linewidth,'markersize',markersize);
set(gca,'box','on','fontname',fontname,'fontsize',fontsize)
xlim([0 100]);
ylim([0 100]);
ylabel('% Targets Viewed')
xlabel('% Distance traveled (normalized)')
title('Search Efficiency with Full System');
% legend([subjectstrings, {'View All'}],'Location','SouthEast')

%% Get Channel Locations
chanlocs = cell(1,numel(subjects));
for i=1:numel(subjects)
    foo = load(sprintf('ALLEEG%d_eyeposcorrected.mat',subjects(i)));
    chanlocs{i} = foo.ALLEEG(1).chanlocs;
end
%% Plot average forward models
bintimes = 100:100:1000; % in ms
binwidth = 100; % in ms
figure;
% figure(255); clf;
% Reconcile different channels being used (use subset common to all subjs)
channies = [chanlocs{:}];
chans = unique({channies.labels});
for i=1:numel(subjects)
    chans = chans(ismember(chans,{chanlocs{i}.labels}));
end
chanlocs_all = chanlocs{1}(ismember({chanlocs{1}.labels},chans));
% chanmatrix = [1:10, 0; 0 0, 11:17, 0, 0; 18:27, 0; 28:38; 39:40,0,41:48; 0,0,0,0,49:51,0,0,0,0; 0,0,0,0,52:54,0,0,0,0];

fwds = zeros(numel(chans), size(R(1).fwdModel,2),numel(subjects));
twts = zeros(size(R(1).v,1), size(R(1).v,2),numel(subjects));
for i=1:numel(subjects)
   fwds(:,:,i) = mean(R(i).fwdModel(ismember({chanlocs{i}.labels},chans),:,:),3);
   twts(:,:,i) = mean(R(i).v,3);
end

% Plot scalp maps
% cmax = max(max(abs(mean(fwds,3))));
nBins = numel(bintimes);
nCols = ceil(nBins/2);
cmax = 10;
for i=1:size(fwds,2)
    subplot(4,nCols,rem(i-1,2)*nCols+floor((i-1)/2)+1,'fontname',fontname,'fontsize',fontsize)
%     topoplot(mean(fwds(:,i,:),3),chanlocs_all,'maplimits',[-cmax cmax],'plotrad',0.5); % don't include mullet
    topoplot(mean(fwds(:,i,:),3),chanlocs_all,'maplimits',[-cmax cmax],'conv','on'); % do include mullet
    title(sprintf('%d-%dms',bintimes(i),bintimes(i)+binwidth));
end
% Annotate plot
axes('Position',[.85 .55 .1 .24],'CLim',[-cmax cmax],'visible','off','fontname',fontname,'fontsize',fontsize);
colorbar('fontname',fontname,'fontsize',fontsize)
% Annotate figure
h = MakeFigureTitle('Forward Models (Mean Across Subjects) (uV)');
set(h,'fontname',fontname,'fontsize',fontsize,'fontweight','normal')

% mode = 'mean';
mode = 'subjects';

subplot(2,1,2,'fontname',fontname,'fontsize',fontsize); hold on;
xbin = bintimes(1:size(twts,2)) + binwidth/2;

switch mode
    case 'mean'
        % Plot temporal weights
        allwts = permute(twts,[2 3 1]);
        errorbar(xbin,mean(allwts,2),std(allwts,[],2),...
            '.-','linewidth',linewidth,'markersize',markersize);        

    case 'subjects'
        % Plot temporal weights
        plot(xbin,permute(twts,[2 3 1]),'.-','linewidth',linewidth,'markersize',markersize);

        % Make legend
        legendstr = cell(1,numel(subjects));
        for i=1:numel(subjects)
            legendstr{i} = sprintf('S%d',subjects(i));
        end
        legend(legendstr, 'Location','SouthWest')

end

% Annotate plot
set(gca,'xgrid','on','box','on')
ylabel('temporal weights (A.U.)')
xlabel('time between saccade offset and bin center (ms)')
hold on
plot(get(gca,'xlim'),[0 0],'k--');
