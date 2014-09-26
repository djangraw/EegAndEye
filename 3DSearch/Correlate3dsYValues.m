function Correlate3dsYValues(R1,R2)

% Correlate3dsYValues(R1,R2)
%
% Created 7/29/13 by DJ.

nSubj = numel(R1);
nPlots = nSubj+2;
nRows = ceil(sqrt(nPlots));
nCols = ceil(nPlots/nRows);
[r,r0,r1,p,p0,p1] = deal(nan(nSubj+1,1));
%% Single-subjects
for i=1:nSubj
    
    [r_tmp,p_tmp] = corrcoef(R1(i).y,R2(i).y); % all trials
    [r0_tmp,p0_tmp] = corrcoef(R1(i).y(R1(i).truth~=1),R2(i).y(R2(i).truth~=1)); % distractors
    [r1_tmp,p1_tmp] = corrcoef(R1(i).y(R1(i).truth==1),R2(i).y(R2(i).truth==1)); % targets

    r(i) = r_tmp(1,2);
    r0(i) = r0_tmp(1,2);
    r1(i) = r1_tmp(1,2);
    p(i) = p_tmp(1,2);
    p0(i) = p0_tmp(1,2);
    p1(i) = p1_tmp(1,2);
    
    subplot(nRows,nCols,i);
    cla; hold on;
    plot(R1(i).y(R1(i).truth~=1),R2(i).y(R2(i).truth~=1),'b.');
    plot(R1(i).y(R1(i).truth==1),R2(i).y(R2(i).truth==1),'r.');
    xlabel('y1');
    ylabel('y2');
    title(sprintf('S%d:\np = %.3g, p0 = %.3g, p1 = %.3g',i,p(i),p0(i),p1(i)));
    
end


%% All subjects 
y1 = [R1.y]';
y2 = [R2.y]';
truth = [R1.truth]';

% calculate
[r_tmp,p_tmp] = corrcoef(y1,y2); % all trials
[r0_tmp,p0_tmp] = corrcoef(y1(truth~=1),y2(truth~=1)); % distractors
[r1_tmp,p1_tmp] = corrcoef(y1(truth==1),y2(truth==1)); % targets

% extract
i = nSubj+1;
r(i) = r_tmp(1,2);
r0(i) = r0_tmp(1,2);
r1(i) = r1_tmp(1,2);
p(i) = p_tmp(1,2);
p0(i) = p0_tmp(1,2);
p1(i) = p1_tmp(1,2);

% Display correlation coefficients
% fprintf('R = %.3g, p = %.3g\n',r(i),p(i));
% fprintf('R0 = %.3g, p0 = %.3g\n',r0(i),p0(i));
% fprintf('R1 = %.3g, p1 = %.3g\n',r1(i),p1(i));

% Plot linear fit
subplot(nRows,nCols,nSubj+1); cla; hold on;
plot(y1(truth~=1),y2(truth~=1),'b.');
plot(y1(truth==1),y2(truth==1),'r.');
xlabel('y1');
ylabel('y2');
title(sprintf('%d subjects:\np = %.3g, p0 = %.3g, p1 = %.3g',nSubj, p(i),p0(i),p1(i)));

% Plot p values, corr coeffs
subplot(nRows,nCols,nSubj+2); cla; hold on;
plot(1:nSubj+1,[r0 r r1],'.-');
plot(find(p0<0.05)-0.2, zeros(1,sum(p0<0.05)), 'b*');
plot(find(p1<0.05)+0.2, zeros(1,sum(p1<0.05)), 'r*');
plot(find(p<0.05), zeros(1,sum(p<0.05)), '*','Color', [0 0.5 0]);
xlabel('Subject')
ylabel('R');
legend('Distractors','Both','Targets');
title('Correlation Coefficients');
