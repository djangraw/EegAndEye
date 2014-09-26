function [contrastFns, group_RF, contrastZ, group_z] = RunTopLevelGlm_Full(R,plus_events,minus_events,baseline_win,multcorrect)

% [contrastFns, group_RF, contrastZ, group_z] =
% RunTopLevelGlm_Full(R,plus_events,minus_events,baseline_win,multcorrect)
%
% Created 4/24/13 by DJ.
% Updated 4/19/14 by DJ - tResponse in cells.

if nargin<4
    baseline_win = [];
end
if nargin<5
    multcorrect = 'none';
end


%% Set up
lastfilename = R(1).filenames{end};
dashes = strfind(lastfilename,'-');
glmType = lastfilename((dashes(3)+1):(length(lastfilename)-4));

%% Set up betas
[contrastFns, contrastVar, contrastZ] = SetUpTopLevelGlm(R,plus_events,minus_events,baseline_win);

for i=1:numel(R)
    figure(100+i); subplot(3,1,1); ylim([-0.5 0.5]);
end

%% Run top-level GLM
[group_RF,group_P] = RunTopLevelGlm_EEG(contrastFns,contrastVar,multcorrect);

%% Plot results
if iscell(R(1).tResponse)
    tResponse = R(1).tResponse{R(1).iLevel};
else
    tResponse = R(1).tResponse;
end
figure;
MakeFigureTitle(sprintf('%s GLM, %s - %s contrast',glmType,plus_events,minus_events));

subplot(2,2,1); hold on;
plot(tResponse, mean(contrastFns,4));
plot(tResponse, mean(mean(contrastFns,4),1),'k','linewidth',2);
xlabel('time (ms)')
ylabel('voltage (uV)')
title('Subject Average')

subplot(2,2,2); hold on;
plot(tResponse, group_RF);
plot(tResponse, mean(group_RF,1),'k','linewidth',2);
xlabel('time (ms)')
ylabel('voltage (uV)')
title('Hierarchical GLM result')

subplot(2,2,3); hold on;
plot(tResponse, mean(contrastFns,4)-group_RF);
plot(tResponse, mean(mean(contrastFns,4)-group_RF,1),'k','linewidth',2);
xlabel('time (ms)')
ylabel('voltage (uV)')
title('Difference')

subplot(2,2,4); hold on;
plot(tResponse, group_P);
plot(tResponse, mean(group_P,1),'k','linewidth',2);
isSignifHi = any(group_P>0.975);
isSignifLo = any(group_P<0.025);
plot(get(gca,'xlim'),[0 0],'k--')
plot(get(gca,'xlim'),[1 1],'k--')
plot(get(gca,'xlim'),[0.025 0.025],'k:')
plot(get(gca,'xlim'),[.975 .975],'k:')
plot(tResponse(isSignifHi),repmat(1.1,1,sum(isSignifHi)),'g.');
plot(tResponse(isSignifLo),repmat(-0.1,1,sum(isSignifLo)),'g.');
ylim([-0.1 1.1]);
xlabel('time (ms)')
ylabel('P value')
title(sprintf('Hierarchical GLM p values (mc correction: %s)',multcorrect))


%% Make z score movie
% Convert to Z scores
group_z = norminv(group_P);

