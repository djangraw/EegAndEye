function PlotRfsWithPermuationStats(b_data,b_rand)

% PlotRfsWithPermuationStats(b_data,b_rand)
%
% Created 8/8/14 by DJ.

[D,Nh,Nr] = size(b_data);
b_rand_all = cat(4,b_rand{:});
pVal = zeros(size(b_data));
% sigDiff = zeros(size(b_data));
for i=1:D
    fprintf('i=%d/%d...\n',i,D);
    for j=1:Nh
        for k=1:Nr   
            if ismember(k,iEvents)
                pVal(i,j,k) = mean((b_data(i,j,k)+b_data(i,j,k+7))>(squeeze(b_rand_all(i,j,k,:)) + squeeze(b_rand_all(i,j,k+7,:)))); 
            else
                pVal(i,j,k) = mean(b_data(i,j,k)>squeeze(b_rand_all(i,j,k,:))); 
            end
%             if mean(b_data(i,j,k)>squeeze(b_rand_all(i,j,k,:)))>0.95 
%                 sigDiff(i,j,k) = 1;
%             elseif mean(b_data(i,j,k)<squeeze(b_rand_all(i,j,k,:)))>0.95
%                 sigDiff(i,j,k) = -1;
%             end
        end
    end
end

sigDiff = pVal>0.95 | pVal<0.05;

%%
chansToPlot = {'FZ';'CZ';'PZ';'OZ'};
event_types = R.regressor_events{R.iLevel};
Cmap = GetSquaresEventColormap(event_types);
chanlocs = R.EEG.chanlocs;
iEvents = 4:7;
% PlotResponseFnsGrid(b_data(:,:,iEvents),event_types(iEvents),tResponse,chanlocs,chansToPlot,Cmap(iEvents,:));
PlotResponseFnsGrid(b_data(:,:,iEvents)+b_data(:,:,iEvents+7),event_types(iEvents),tResponse,chanlocs,chansToPlot,Cmap(iEvents,:));

iChan = nan(size(chansToPlot));
for i=1:numel(chansToPlot)
    iChan(i) = find(strcmp( chansToPlot{i},{chanlocs.labels}));
end
for i=1:numel(chansToPlot)
    subplot(numel(chansToPlot),1,i);
    for j=1:numel(iEvents)
        starTimes = tResponse(sigDiff(iChan(i),:,iEvents(j))>0);
        plot(starTimes,repmat(1.5+.1*j,size(starTimes)),'*','Color',Cmap(iEvents(j),:));
    end
end
    
