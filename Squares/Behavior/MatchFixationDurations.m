function [isKeeper, class, squarefixtime] = MatchFixationDurations(y,classes,tLong,reject_behavior_rules,method,histedges)

% [isKeeper, class, squarefixtime] = MatchFixationDurations(y,classes,tLong,reject_behavior_rules,method,histedges)
% 
% INPUTS:
% -y is a behavioral dataset from loadBehaviorData.
% -classes is a cell array of strings, to be used as input to UseEventRule.
% -tLong is a scalar indicating the time (in ms) at which the distance
% between events will cause their response functions to overlap.
% (If the response function goes from tStart to tEnd relative to the event,
% tLong = tEnd-tStart).
% -reject_behavior_rules is a cell array of strings to be used as input to
% RejectBehavior.
% -histedges is an array of times (in ms) that you want to use as the edges
% of your histogram bins in the output plots.
%
% OUTPUTS:
% -isKeeper is an Nx5 matrix of booleans indicating whether the saccade to
% each square (row=trial, column=sqnum) should be kept in the matched set.
% -class is an Nx5 cell array of strings indicating the class of each
% square (only those in 'classes' are included).
% -squarefixtime is an Nx5 matrix of times in ms that the subject fixated
% on the square.
%
% Created 4/26/13 by DJ.

% Handle inputs
if nargin<3 || isempty(tLong)
    tLong = 500; % in ms
end
if nargin<4
    reject_behavior_rules = ''; % no behavior is too bad!
end
if nargin<5 || isempty(method)
    method = 'KSTest';
end
if nargin<6 || isempty(histedges)
    % Declare histogram bins
    histedges = 0:100:1000; % in ms
end


% --- STEP 0: Get all fixation durations

% Declare constants
nClasses = length(classes);
nBins = numel(histedges)-1;
colors = {'r' 'g' 'b' 'c' 'm' 'y' 'r--' 'g--' 'b--' 'c--' 'm--' 'y--'};

% set up
nSessions = numel(y);
nTrialsPerSession = numel(y(1).trial.start_time);
nTrials = nTrialsPerSession*nSessions;

class = repmat({''},nTrials,5);
squarefixtime = nan(nTrials,5);
isKeeper = true(nTrials,5);

disp('Classifying events...')

for k=1:nSessions
    s = y(k).saccade;
    % Get duration of fixating on this square
    iChangeSquare = find(diff(s.squarenum)>0 | diff(s.squarenum)<0)+1; % saccades when subject moved away from current square    
    squaredur = [s.start_time(iChangeSquare(2:end)) - s.end_time(iChangeSquare(1:end-1)); inf];
    isNotEnd = (s.squarenum(iChangeSquare)>0 & s.squarenum(iChangeSquare)<6);
    iChangeSquare = iChangeSquare(isNotEnd);
    squaredur = squaredur(isNotEnd);
    % get trial and squarenum
    trial = s.trialnum(iChangeSquare) + (k-1)*nTrialsPerSession;
    sqnum = s.squarenum(iChangeSquare);
    for j=1:numel(squaredur)
        squarefixtime(trial(j),sqnum(j)) = squaredur(j);
    end
    
end

for i=1:nClasses
    [~,~,~,iEvents] = UseEventRule(y,classes{i});
    for k=1:nSessions     
        s = y(k).saccade;
        trial = s.trialnum(iEvents{k}) + (k-1)*nTrialsPerSession;
        sqnum = s.squarenum(iEvents{k});
        for j=1:numel(iEvents{k})
            if sqnum(j)>0 & sqnum(j)<6
                class{trial(j),sqnum(j)} = classes{i};
            end
        end
        
    end
end
% Reject bad trials
for k=1:nSessions
    isBadTrial = RejectBehavior(y(k),reject_behavior_rules);
    badTrials = find(isBadTrial) + (k-1)*nTrialsPerSession;
    isKeeper(badTrials,:) = false;
end

% If no method, just stop here
if strcmp(method,'none')
    return;
end

% Group into chunks of events with overlapping responses (using tLong cutoff)
chunk = ones(nTrials,5);
isLong = squarefixtime>tLong;
chunk(:,2:5) = 1 + cumsum(isLong(:,1:4),2);
for i=2:nTrials
    chunk(i,:) = chunk(i,:) + chunk(i-1,end);
end

% get Before Histogram
subplot(2,2,1);
is1 = (strcmp(classes{1},class) & isKeeper);
is2 = (strcmp(classes{2},class) & isKeeper);
[~,p_0] = kstest2(squarefixtime(is1),squarefixtime(is2)); 
classhist = GetClassHist(classes,class,squarefixtime,isKeeper,histedges);
nPerClass_before = sum(classhist,2);
title(sprintf('Before matching: p = %g',p_0))

go = true;
lockClass = false(nClasses,1);
z_after = norminv(p_0); % initialize
n = [];

while go
    switch method
        case 'BinSubtract'
            % find smallest class
            nPerClass = sum(classhist,2);
            nInSmallestClass = min(nPerClass(~lockClass));
            iSmallestClass = find(nPerClass==nInSmallestClass & ~lockClass,1);
            fprintf('%s is smallest class (%d examples)\n',classes{iSmallestClass},nInSmallestClass);        
            % if this class has any bins where it's bigger than another class,
            % remove an offending trial
            binisok = false(1,nBins);
            for j=1:nBins
                if any(classhist(:,j)<classhist(iSmallestClass,j))
                    squaresInBin = find(strcmp(class,classes{iSmallestClass}) & squarefixtime>=histedges(j) & squarefixtime<histedges(j+1) & isKeeper);
                    chunksInBin = chunk(squaresInBin);
                    nToRemove = max(classhist(iSmallestClass,j)-classhist(:,j));
                    fprintf('removing %d chunks from bin %d\n',nToRemove,j);
                    % Just pick 1st N
                    isKeeper(ismember(chunk,chunksInBin(1:nToRemove))) = false;          
        %             isKeeper(trialsInBin(1:nToRemove),:) = false;
                    break % exit loop to recalculate histogram
                else
                    binisok(j) = true;
                end
            end   
            if all(binisok)
                disp('class done!')
                lockClass(iSmallestClass) = true;
            end    
            classhist = GetClassHist(classes,class,squarefixtime,isKeeper,histedges);
            if all(lockClass)
                go = false;
            end
            
        case 'KSTest'
            if numel(classes)>2
                error('KSTest option only works for 2 classes!');
            end
            is1 = (strcmp(classes{1},class) & isKeeper);
            is2 = (strcmp(classes{2},class) & isKeeper);
            n = cat(2,n,[sum(is1(:)); sum(is2(:))]);

            [~,p_before] = kstest2(squarefixtime(is1),squarefixtime(is2));            
            % Find chunk that will improve p_before the most
            chunksInPlay = unique(chunk(is1 | is2));
            fprintf('%d chunks in play...',numel(chunksInPlay));
            p_new = zeros(1,numel(chunksInPlay));
            for j=1:numel(chunksInPlay)
                [~,p_new(j)] = kstest2(squarefixtime(is1 & chunk~=chunksInPlay(j)),squarefixtime(is2 & chunk~=chunksInPlay(j))); 
            end
            [p_after,iRemove] = max(p_new);
            if p_after<=p_before
                break
            else
                z_after = cat(2,z_after, norminv(p_after));
                fprintf('removed 1 to make Z = %g\n',z_after(end));
                isKeeper(chunk==chunksInPlay(iRemove)) = false;
%                 % plot histogram
%                 classhist = GetClassHist(classes,class,squarefixtime,isKeeper,histedges);
            end

    end
end


% get After Histogram
disp('SUCCESS!');
subplot(2,2,2)
classhist = GetClassHist(classes,class,squarefixtime,isKeeper,histedges);
title(sprintf('After matching: p = %g',p_after))
set(subplot(2,2,2),'ylim',get(subplot(2,2,1),'ylim'));
nPerClass_after = sum(classhist,2);
disp('---BEFORE Matching:---')
fprintf('KSTest: p = %g\n',p_0);
for i=1:numel(classes)
    fprintf('%s: %d examples\n',classes{i},nPerClass_before(i));
end
disp('---AFTER Matching:---')
fprintf('KSTest: p = %g\n',p_after);
for i=1:numel(classes)
    fprintf('%s: %d examples\n',classes{i},nPerClass_after(i));
end
disp('---------------------');

% Plot z timecourse
subplot(2,1,2);
h_ax = plotyy(0:length(z_after)-1, n, 0:length(z_after)-1, z_after);
xlabel('iterations');
ylabel(h_ax(1),'nInClass');
legend(h_ax(1),classes,'Location','North');
ylabel(h_ax(2),'Z');
% ylabel('KSTest z score');





% HELPER FUNCTION
function classhist = GetClassHist(classes,class,squarefixtime,isKeeper,histedges)

classhist = nan(numel(classes),numel(histedges));
for i=1:numel(classes)
    classhist(i,:) = histc(squarefixtime(strcmp(classes{i},class) & isKeeper),histedges);
end
% Get histograms of fixation times for each class
diffedges = diff(histedges);
histcenters = histedges + [diffedges diffedges(end)]/2;
plot(histcenters,classhist','linewidth',2);
grid on
% bar(histedges,classhist','histc');

% Annotate plot
xlabel('Duration of fixation on this square (ms)');
ylabel('# fixations');
legend(classes);
% % More informative legend
% nPerClass = sum(classhist,2);
% legendtext = cell(1,numel(classes));
% for i=1:numel(classes)
%     legendtext{i} = sprintf('%s (n=%d)',classes{i},nPerClass(i));
% end
% legend(legendtext);


