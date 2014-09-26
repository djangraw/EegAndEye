for iSubj=1:numel(subjects)
    fprintf('Subject %d...\n',iSubj);
    R(iSubj) = load(sprintf('%s-%d-Type-v3pt6-RampUp-RidgeTrace',experiment,subjects(iSubj)));    
end
%% Get covariance matrices
XTX_all = [];
for i=1:numel(R)
    XTX_all = cat(3,XTX_all,R(i).XTX);
end

%% Plot covariance matrices
% event_list = Rnew.regressor_events{1}(2:end);
tResponse = 0:10:1000;
T = length(tResponse);

subplot(1,2,1);
imagesc(mean(XTX_all,3));
set(gca,'xtick',(1:T:size(XTX_all,1)),'xticklabel',event_list);
set(gca,'ytick',(1:T:size(XTX_all,1)),'yticklabel',event_list);
rotateticklabel(gca,90);
colorbar;
title(sprintf('Mean covariance matrix across %d %s subjects',numel(subjects),experiment))
set(gca,'clim',[0 0.03])
grid on;

subplot(1,2,2);
imagesc(var(XTX_all,0,3));
set(gca,'xtick',(1:T:size(XTX_all,1)),'xticklabel',event_list);
set(gca,'ytick',(1:T:size(XTX_all,1)),'yticklabel',event_list);
rotateticklabel(gca,90);
colorbar;
set(gca,'clim',[0 0.001])
title(sprintf('Variance in covariance matrices across %d %s subjects',numel(subjects),experiment))
grid on;



