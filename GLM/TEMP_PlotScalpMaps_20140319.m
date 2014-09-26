
% Set up
prefix = 'sf3';
eval(sprintf('R = R_%s_type;',prefix));
R = UpdateGlmResultsFormat(R);

switch prefix
    case 'sf3'
        event_pairs = {'pD_{1/3}','pD_{0/3}';...
            'pD_{2/3}','pD_{0/3}';
            'pD^*_{-/3}','pD_{0/3}'};
        
%         event_pairs = {'pT_{0/3}','pD_{0/3}';...
%             'pT_{1/3}','pD_{0/3}';...
%             'pT^*_{2/3}','pD_{0/3}'};
%         event_pairs = {'pD_{1/3}','pD_{2/3}';...
%             'pD_{0/3}','pD_{1/3}';...
%             'pT_{0/3}','pD_{0/3}';...
%             'pT_{1/3}','pT_{0/3}';...
%             'pT^*_{2/3}','pT_{1/3}'};
    case 'sf'
        event_pairs = {'pD_{1/2}','pD_{0/2}';...            
            'pD^*_{-/2}','pD_{0/2}';...
            'pD_{+/2}','pD_{0/2}'};
%         event_pairs = {'pT_{0/2}','pD_{0/2}';...            
%             'pT^*_{1/2}','pD_{0/2}';...
%             'pT_{+/2}','pD_{0/2}'};
%         event_pairs = {'pD_{0/2}','pD_{1/2}';...
%             'pT_{0/2}','pD_{0/2}';...
%             'pT^*_{1/2}','pT_{0/2}'};
    case 'sq'
        event_pairs = {'aD_{1/2}','aD_{0/2}';...
            'aD^*_{-/2}','aD_{0/2}';...
            'aD_{+/2}','aD_{0/2}'};
%         event_pairs = {'aT_{0/2}','aD_{0/2}';...
%             'aT^*_{1/2}','aD_{0/2}';...
%             'aT_{+/2}','aD_{0/2}'};
%         event_pairs = {'aD_{1/2}','aD_{0/2}';...
%             'aT_{0/2}','aD_{0/2}';...
%             'aT^*_{1/2}','aT_{0/2}'};
end

% Calculate contrast
foo = cell(1,numel(R));
for i=1:numel(R)
    for j=1:size(event_pairs,1)
        i1 = strcmp(event_pairs{j,1},R(i).regressor_events{3});
        i2 = strcmp(event_pairs{j,2},R(i).regressor_events{3});
        foo{j,i} = R(i).responseFns{3}(:,:,i1) - R(i).responseFns{3}(:,:,i2);
    end
end

% Compile
data = [];
for j=1:size(event_pairs,1)
    data = cat(3,data,mean(cat(4,foo{j,:}),4));
end

%% Plot
% Set up
tBinCenters = 37.5:75:750;
tBinWidth = 75;
clim = [-1.2 1.2]; % make empty for automatic

legendstr = cell(1,size(event_pairs,1));
for j=1:size(event_pairs,1)
    legendstr{j} = sprintf('%s - %s',event_pairs{j,1}, event_pairs{j,2});
end

% Plot scalp maps
% figure;
clf;
[sm_all,sm] = GetScalpMaps(data,R(1).tResponse{end},tBinCenters,tBinWidth);
PlotScalpMaps(sm_all,R(1).EEG.chanlocs,clim,tBinCenters-tBinWidth/2,legendstr);
% Annotate figure
MakeFigureTitle(sprintf('Mean across %d %s Subjects',numel(R),prefix));
