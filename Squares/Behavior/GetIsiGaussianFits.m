function obj = GetIsiGaussianFits(results,class)

% Created 2/5/13 by DJ.

% Declare histogram bins
xIsi = 0:5:500;%10:20:1000;



% Get event indices 
[~,~,~,iEvents] = UseEventRule(results,class);
% Get event info
nSessions = numel(results);
[dist,fixdurA] = deal(cell(1,nSessions));
for k=1:nSessions
    s = results(k).saccade;
    dist{k} = sqrt((s.end_position(iEvents{k},1)-s.start_position(iEvents{k},1)).^2 + ...
        (s.end_position(iEvents{k},2)-s.start_position(iEvents{k},2)).^2);

    fixdur = s.start_time(2:end) - s.end_time(1:end-1);
    fixdurA{k} = fixdur(iEvents{k}); % duration of fixation after this saccade
end
distance = cat(1,dist{:});
fixdur_after = cat(1,fixdurA{:});

%%
% Fit distance distribution to mixture of gaussians (MoG)
X = fixdur_after(fixdur_after<xIsi(end)); %cut off at 1s
% X = fixdur_after(fixdur_after<xIsi(end) & distance>100); %cut off at 1s
nGaussians = 2;
obj = gmdistribution.fit(X,nGaussians,'sharedcov',false,'replicates',10);

% Get gaussian distributions
amp = obj.PComponents;
mu = [obj.mu(:)];
sigma = sqrt([obj.Sigma(:)]);
%% Plot results
cla; hold on;
nIsi = hist(X,xIsi);
bar(xIsi,nIsi,'facecolor','k');
foo = zeros(nGaussians,length(xIsi));
for i=1:nGaussians
    fprintf('Gaussian #%d: (A,mu,sigma)  = (%.2g, %.0f, %.0f)\n',i,amp(i),mu(i),sigma(i));
    foo(i,:) = amp(i)*normpdf(xIsi,mu(i),sigma(i));       
end
normfactor = sum(nIsi)/sum(foo(:));
plot(xIsi,foo*normfactor);    

% plot(xIsi,sum(foo,1)*normfactor,'k');
foo = obj.pdf(xIsi');
plot(xIsi,foo/sum(foo)*sum(nIsi),'m--');
xlabel('fixation duration (ms)')
ylabel('# fixations')
legend(sprintf('hist of ISI<%dms',xIsi(end)),'1st gaussian','2nd gaussian','fit','Location','NorthWest');