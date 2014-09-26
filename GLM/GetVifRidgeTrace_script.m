
%% If running RunEegGlm, get Scrop and then run this cell.

p = size(Scrop,2);
mnS = zeros(1,p);
for i=1:p
    Scrop(:,i) = (Scrop(:,i) - mnS(i)) / sqrt(sum((Scrop(:,i) - mnS(i)).^2));
end
    
XTX = Scrop'*Scrop;

%% Sweep k/lambda values and find max VIF for each
% ks = 0:0.001:0.03; % for SF experiments
ks = 0:.002:.02;%0:1e-8:1e-7;%0:1e-6:1e-6; % for SQ
VIFmax = zeros(size(ks));
for i=1:numel(ks)
    fprintf('%d/%d...\n',i,numel(ks));
    k = ks(i);
    VIFmax(i) = max(diag((XTX+k*eye(p))^(-1)*XTX*(XTX+k*eye(p))^(-1)));
end
    
%% Plot results
% Plot VIF vs. k
cla;
plot(ks,VIFmax)
xlabel('ridge bias parameter')
ylabel('Variance Inflation Factor')
title('VIF vs. k for sq subject 3')

% Find lowest k for which VIF<10.
hold on;
plot(get(gca,'xlim'),[10 10],'k--');
iBest = find(VIFmax<=10,1);
plot([ks(iBest) ks(iBest)],[0 VIFmax(iBest)],'k--');
ylim([0 50]);