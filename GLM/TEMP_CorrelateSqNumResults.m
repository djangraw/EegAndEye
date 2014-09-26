% TEMP_CorrelateSqNumResults.m
%
% Created 1/23/13 by DJ.

% First run GetGroupSvdResults_script on Type-SqNum data
tWin = [0 500]; % time window in ms
if iscell(S(1).tResponse)
    tResponse = S(1).tResponse{S(1).iLevel};
else
    tResponse = S(1).tResponse;
end
oktime = tResponse>tWin(1) & tResponse<tWin(2);
colors = 'rgbcm';
cla; hold on
disp('-----')
for i=2:5
    
    scatter(avgTimecourses(1,oktime,1),avgTimecourses(i,oktime,1),[colors(i) '.'])
end
for i=2:5
%     b = regress(avgTimecourses(i,oktime,1)',[avgTimecourses(1,oktime,1)',ones(sum(oktime),1)]); % with offset
%     b = [regress(avgTimecourses(i,oktime,1)',avgTimecourses(1,oktime,1)') 0]; % without offset
    b = [-1/2^(i-1), 0]; % assume 1/2 rolloff
    plot(avgTimecourses(1,oktime,1),avgTimecourses(1,oktime,1)*b(1)+b(2),colors(i));
%     plot(avgTimecourses(1,oktime,1),avgTimecourses(1,oktime,1)*-1/2^(i-1),colors(i));
    fprintf('SN%d = SN1 * %.2f + %.2f\n',i,b(1),b(2));

end

xlabel('SqNum1 voltage')
ylabel('SqNumX voltage')
legend('SqNum2','SqNum3','SqNum4','SqNum5');
title(sprintf('%d-%d ms',tWin(1),tWin(2)));

%% Linear regression
cla; hold on;
plot(tResponse(oktime),avgTimecourses(1,oktime,1),colors(1));
disp('-----')
for i=2:5
%     b = regress(avgTimecourses(1,oktime,1)',[avgTimecourses(i,oktime,1)',ones(sum(oktime),1)]);
    b = [regress(avgTimecourses(1,oktime,1)',avgTimecourses(i,oktime,1)') 0];
    plot(tResponse(oktime),avgTimecourses(i,oktime,1)*b(1)+b(2),colors(i));
    fprintf('SN1 = SN%d * %.2f + %.2f\n',i,b(1),b(2));
end
xlabel('time (ms)')
ylabel('voltage (V)')
legend('SqNum1','SqNum2','SqNum3','SqNum4','SqNum5');
title(sprintf('%d-%d ms',tWin(1),tWin(2)));

%%
cla; hold on;
mults = [1 -1/2 -1/4 -1/8 -1/16];

for i=1:5
    plot(tResponse(oktime),avgTimecourses(i,oktime,1),colors(i));
    plot(tResponse(oktime),avgTimecourses(1,oktime,1) * mults(i),[colors(i) '--']);
end
xlabel('time (ms)')
ylabel('voltage (V)')
legend('SqNum1',sprintf('SN1*%.2f',mults(1)),'SqNum2',sprintf('SN1*%.2f',mults(2)),'SqNum3',sprintf('SN1*%.2f',mults(3)),'SqNum4',sprintf('SN1*%.2f',mults(4)),'SqNum5',sprintf('SN1*%.2f',mults(5)));
title(sprintf('%d-%d ms',tWin(1),tWin(2)));
