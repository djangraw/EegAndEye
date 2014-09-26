% Check For Offset
%
% Load up the GLM results across subjects and plot each one's "average ERP"
% So they can be compared to each other side by side.
%
% Created 3/5/14 by DJ.

%% Set up

% prefix = 'sq';
% subjects = [9:11 13:15 17:19];
% % glmTypes = {'LREvents-SqNum-Type-v2pt4'};
% glmTypes = {'LREvents-Type-SqNum-v2pt4'};

prefix = 'sf';
subjects = [1:10 12];
% glmTypes = {'SqNum-Type-v2pt4'};
glmTypes = {'SqNum-Type-v3pt0'};

% prefix = 'sf3';
% subjects = 1:12;
% glmTypes = {'SqNum-Type-v3pt0'};

% bandnames = {'', 'theta_v2' 'lower1alpha_v2' 'lower2alpha_v2' 'upperalpha_v2','delta_v2' 'beta_v2' 'lowergamma_v2'};
% bandnames = {'', 'Theta_v2-' 'Lower1Alpha_v2-' 'Lower2Alpha_v2-' 'UpperAlpha_v2-','Delta_v2-' 'Beta_v2-' 'LowerGamma_v2-'};
bandnames = {''};

%% Load

% Load
clear R
bandname = bandnames{1};
glmType = glmTypes{1};
for j=1:numel(subjects)
    fprintf('Loading subject %d/%d...\n',j,numel(subjects));
    filename = sprintf('%s-%d-GLMresults-%s%s',prefix,subjects(j),bandname,glmType);
    foo = load(filename);
    if j==1
        R(j) = foo;
    else
        R(j) = orderfields(foo,R(1));
    end
end
disp('Done!')
    
%% Plot Average ERPs
if cell(R(1).tResponse)
    tResponse = R(1).tResponse{R(1).iLevel};
else
    tResponse = R(1).tResponse;
end
data = [];
legendstr = repmat({''},1,numel(subjects));
for j = 1:numel(subjects)
    data(:,:,j) = mean(R(j).responseFns{1}(:,:,1:5),3); % mean across squares event types
    legendstr{j} = sprintf('Subject %d',subjects(j));
end

Cmap = distinguishable_colors(numel(subjects),{'w','k'}); % don't allow black or white

figure(111); MakeFigureTitle('Average ERP (batch 1)');
PlotResponseFnsGrid(data(:,:,1:6),legendstr(1:6),tResponse,R(iSubj).EEG.chanlocs,{'FZ'; 'CZ'; 'PZ'},Cmap(1:6,:));
figure(112); MakeFigureTitle('Average ERP (batch 2)');
PlotResponseFnsGrid(data(:,:,7:end),legendstr(7:end),tResponse,R(iSubj).EEG.chanlocs,{'FZ'; 'CZ'; 'PZ'},Cmap(7:end,:));



%% Plot Other Events

nEventTypes = 6;
for i=1:nEventTypes;
    eventdata = [];
    legendstr = repmat({''},1,6);%numel(subjects));
    for j = 1:numel(subjects)
        eventdata(:,:,j) = data(:,:,j) + R(j).responseFns{2}(:,:,i); % mean across squares event types
        legendstr{j} = sprintf('Subject %d',subjects(j));
    end

    Cmap = distinguishable_colors(numel(subjects),{'w','k'}); % don't allow black or white

    figure(120 + i); MakeFigureTitle(R(1).regressor_events{2}{i});
    PlotResponseFnsGrid(eventdata(:,:,1:6),legendstr(1:6),tResponse,R(iSubj).EEG.chanlocs,{'FZ'; 'CZ'; 'PZ'},Cmap(1:6,:));
    figure(130 + i); MakeFigureTitle(R(1).regressor_events{2}{i});
    PlotResponseFnsGrid(eventdata(:,:,7:end),legendstr(7:end),tResponse,R(iSubj).EEG.chanlocs,{'FZ'; 'CZ'; 'PZ'},Cmap(7:end,:));
    
end
linkaxes(GetSubplots( [120+(1:nEventTypes), 130+(1:nEventTypes)]));

%% Plot Other Events
% Set up
eventdata = [];
legendstr = R(1).regressor_events{end}(1:nEventTypes);
Cmap = distinguishable_colors(nEventTypes,{'w','k'}); % don't allow black or white
% Get
for i=1:nEventTypes;
    
    for j = 1:numel(subjects)
        eventdata(:,:,i,j) = data(:,:,j) + R(j).responseFns{2}(:,:,i); % mean across squares event types
    end    
end
% Plot
figure(140); MakeFigureTitle([prefix '-' glmTypes{end},bandnames{end}]);
PlotResponseFnsGrid(eventdata,legendstr,tResponse,R(1).EEG.chanlocs,{'FZ'; 'CZ'; 'PZ'},Cmap);    

