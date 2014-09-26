% PlotTrialTypeExpectations
%
% -Considering all the possible permutations of red, green, and blue 
% squares in the 5-square experiment (variable 'perms'), this program plots
% the probability that each trial will have z targets among the first x 
% squares seen (variable 'fracTrials').
% -It also plots the probability that each trial will be a complete trial,
% given that the subject has seen z targets among the first x squares
% (variable 'compProb').
%
% Created 1/12/12 by DJ.

% Get all possible permutations of the 3 colors
perms = zeros(3^5,5);
iRow = 0;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for m=1:3
                    iRow = iRow+1;
                    perms(iRow,:) = [i j k l m];
                end
            end
        end
    end
end

% Find fraction of trials that meet certain criteria
isTarg = perms==1; % 1 for targets, 0 for all others
isComp = sum(isTarg,2)>1; % is each trial a 'complete' trial?
fracTrials = nan(6,6); % fraction of trials in each condition
compProb = nan(6,6); % probability that trials in this condition will be complete trials
fracType = nan(6,5);
for i=0:5
    if i>0
        fracType(i,1) = mean(~isTarg(:,i) & sum(isTarg(:,1:i-1),2)==0); % pre-integration distractor
        fracType(i,2) = mean(isTarg(:,i) & sum(isTarg(:,1:i-1),2)==0); % integration target
        fracType(i,3) = mean(~isTarg(:,i) & sum(isTarg(:,1:i-1),2)==1); % post-integration distractor
        fracType(i,4) = mean(isTarg(:,i) & sum(isTarg(:,1:i-1),2)==1); % completion target
        fracType(i,5) = mean(sum(isTarg(:,1:i-1),2)>1); % irrelevant square (target or distractor)
    end
    for j=0:5
        fracTrials(i+1,j+1) = sum(sum(isTarg(:,1:i),2)==j)/size(perms,1);
        compProb(i+1,j+1) = mean(isComp(sum(isTarg(:,1:i),2)==j));        
    end
end
clear i j k l m iRow isTarg isComp

fracType(6,:) = mean(fracType(1:5,:),1);

% Plot results
figure(1)
plot(0:5,fracTrials,'.-')
set(gca,'xtick',0:5)
xlabel('squares seen so far')
ylabel('fraction of trials')
legend('0T so far', '1T so far', '2T so far', '3T so far', '4T so far', '5T so far');

figure(2)
plot(0:5,compProb(:,1:2),'.-')
set(gca,'xtick',0:5)
xlabel('squares seen so far')
ylabel('fraction of trials that are complete')
legend('0T so far', '1T so far');

figure(3)
plot(1:6,fracType,'.-')
set(gca,'xtick',1:6,'xticklabel',{'1','2','3','4','5','all'})
xlabel('squares seen so far')
ylabel('fraction of squares')
legend('D-0T', 'Integ', 'D-1T', 'Compl', 'Irrel');
