% function CompareSquareResponses()
% Created 4/19/14 by DJ.

prefixes = {'sq','sf','sf3'};
iLevel = 3;
events = {'SqNum2','sf-SqNum2','sf-SqNum2'};
tBinCenters = 25:50:500;
tBinWidth = 50;
clim = [-2.5 2.5]; % make empty for automatic

sm = zeros(R_sf_sqnum(1).EEG.nbchan,length(tBinCenters),numel(prefixes));
for iExp = 1:numel(prefixes)
    prefix = prefixes{iExp};

    % Set up
    eval(sprintf('R = R_%s_sqnum;',prefix));
    R = UpdateGlmResultsFormat(R);

    
    % Calculate contrast
    foo = cell(1,numel(R));
    for i=1:numel(R)
        for j=1:size(events,1)            
            i1 = strcmp(events{j,iExp},R(i).regressor_events{iLevel});
%             i2 = strcmp(event_pairs{j,2},R(i).regressor_events{3});
            foo{j,i} = R(i).responseFns{iLevel}(:,:,i1);
%             foo{j,i} = R(i).responseFns{3}(:,:,i1) - R(i).responseFns{3}(:,:,i2);
        end
    end

    % Compile
    data = [];
    for j=1:size(events,1)
        data = cat(3,data,mean(cat(4,foo{j,:}),4));
    end
    [~,sm(:,:,iExp)] = GetScalpMaps(data,R(1).tResponse{end},tBinCenters,tBinWidth);
    
end

%% Plot   
    % Set up
    legendstr = events;

    % Plot scalp maps
    figure(351);
    clf;
%     [sm_all,sm] = GetScalpMaps(data,R(1).tResponse{end},tBinCenters,tBinWidth);
    PlotScalpMaps(sm,R(1).EEG.chanlocs,clim,tBinCenters-tBinWidth/2,legendstr);
    % Annotate figure
%     MakeFigureTitle(sprintf('Mean across %d %s Subjects',numel(R),prefix));