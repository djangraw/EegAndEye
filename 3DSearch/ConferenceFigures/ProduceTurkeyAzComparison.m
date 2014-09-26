% ProduceTurkeyAzComparison
% Created 7/7/11 by DJ

% set up
figure(74)
clf;
% plot
A = [.67 .64 .59; .76 .69 .65];
plot(A)
hold on
% annotate plot
xlim([0.5 2.5]);
set(gca,'xtick',[1 2],'xticklabels',{'Standard','Jittered'})
xlabel('Logistic Regression Type')
ylabel('Maximum Az between -0.5 and 1.0 s')
ylim([0.3 1])
title('Logistic Regression Comparison')
legend('S2','S14','S15')
% plot extra lines
plot(get(gca,'XLim'), [0.5 0.5], 'k--')
plot(get(gca,'XLim'), [0.75 0.75], 'k:')