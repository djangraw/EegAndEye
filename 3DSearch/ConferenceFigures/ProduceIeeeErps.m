% ProduceIeeeErps
%
% Plots the ERPs for Fz, Cz, and Pz for the subjects indicated in the code.
% This includes standard error bars for targets and distractors using the
% JackKnife function, as well as a line for the difference.
% This plot was Figure 2 in DJ's 2011 IEEE/EMBS NE conference submission.
%
% Created 1/14/11 by DJ.
% Updated 1/18/11 by DJ - comments.


% Set up
subjects = [6 2 7];
% subjects = [8 10 11 12];
% channels = [6 15 25]; % Fz Cz Pz, for 'noduds' datasets
channels = [6 16 26]; % Fz Cz Pz, for new 'noduds' datasets

% Load Data
LoadAllEpochs;

% Create and size figure
figure;
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 600 750]);

% Get and plot ERPs
for i = 1:numel(channels)
    % Amass data
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
    
    % Plot Channel ERPs
    subplot(3,1,i); 
    [~, ~, sigTimes] = PlotChannelErp(dataTarg,dataDis,ALLEEG(1).times,ALLEEG(1).chanlocs(channels(i)).labels);
    
    % Reannotate Plot (override PlotChannelErp defaults)
    title('')
    ylabel(sprintf('Channel %s (uV)',ALLEEG(1).chanlocs(channels(i)).labels))
    if i<numel(channels)
        xlabel('');
    end
%     xlim('auto')

end
clear channels dataTarg dataDis i j