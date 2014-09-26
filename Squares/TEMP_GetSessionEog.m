% TEMP_GetSessionEog.m
%
% Created 6/28/12 by DJ for one-time use.

%% Get bahavior data
y = loadBehaviorData(13);
for i=1:10
    EEG13(i) = pop_loadset(sprintf('sq-13-%d-filtered-1000Hz.set',i));
    for j=1:EEG13(i).nbchan
        EEG13(i).data(j,:) = EEG13(i).data(j,:) - mean(EEG13(i).data(j,:));
    end
end

%% Get eog components for each session
for i=1:10
    [component{i}, offset{i}] = FindEogComponent_Squares(EEG13(i),y(i));
end

%% Plot eog components
figure(557);
for i=1:10
    subplot(2,10,i)
    topoplot(component{i},EEG13(i).chanlocs,'electrodes','on');
    colorbar
    title(sprintf('sq-13-%d HEOG component',i))
    subplot(2,10,10+i)
    topoplot(offset{i},EEG13(i).chanlocs,'electrodes','on');
    colorbar
    title(sprintf('sq-13-%d HEOG offset',i))
end
    
    
    
    
    