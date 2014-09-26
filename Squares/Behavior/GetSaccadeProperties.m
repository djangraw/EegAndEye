function GetSaccadeProperties(results,class,reject_rule)

% Plot histograms of saccade distance, duration, and time to previous or
% next saccade.
%
% GetSaccadeProperties(results,class,reject_rule)
%
% INPUTS:
% -results is a cell array in which each cell contains a squares dataset 
% imported by import_squares_data.m.
% -class is a number (see GetSquaresConstants.m) or a string see 
% UseEventRule.m for details) indicating the saccade type to analyze.
% -reject_rule is a string or cell array of strings indicating what rule(s)
% you'd like to use to reject behavior (see RejectBehavior.m for options).
%
% Created 1/5/12 by DJ.
% Updated 5/10/12 by DJ - comments
% Updated 7/16/12 by DJ - fixed colors, mean is now thich black line
% Updated 8/2/12 by DJ - added string support for class input, reject_rule 
%  input

if nargin<3
    reject_rule = '';
end

% Declare constants
nSubjects = length(results);
colors = {'r' 'g' 'b' 'c' 'm' 'y' 'r--' 'g--' 'b--' 'c--' 'm--' 'y--'};
% Declare histogram bins
xDistance = 5:10:300;
xDuration = 2:4:100;
xIsi = 10:20:1000;

% Set up loop
subjects = cell(1,nSubjects);
nDistance = zeros(nSubjects,length(xDistance));
nDuration = zeros(nSubjects,length(xDuration));
nIsiBefore = zeros(nSubjects,length(xIsi));
nIsiAfter = zeros(nSubjects,length(xIsi));
clf; % clear current figure

% Main loop
for i=1:nSubjects
    % Find relevant saccades
%     nSessions = length(results{i});        
    s = combine_structs([results{i}.saccade],[],'columns');
    if isnumeric(class)
        iOk = find(s.class==class);
    else
        [~,~,~,iGood] = UseEventRule(results{i},class);
        for j=1:numel(results{i})
            [~,isBad] = RejectBehavior(results{i}(j),reject_rule);
            iEvents{j} = setdiff(iGood{j},find(isBad));
        end
        iOk = iEvents{1};
        offset = 0;
        for j=2:numel(iEvents)
            offset = offset + numel(results{i}(j-1).saccade.end_time);
            iOk = [iOk; iEvents{j}+offset];
        end
    end
    % Get saccade info
    distance = sqrt((s.end_position(iOk,1)-s.start_position(iOk,1)).^2 + ...
        (s.end_position(iOk,2)-s.start_position(iOk,2)).^2);
    duration = s.end_time(iOk) - s.start_time(iOk);
    isi = diff(s.start_time);
    isi_before = isi(iOk-1);
    isi_after = isi(iOk);
    nDistance(i,:) = hist(distance,xDistance);
    nDuration(i,:) = hist(duration,xDuration);
    nIsiBefore(i,:) = hist(isi_before,xIsi);
    nIsiAfter(i,:) = hist(isi_after,xIsi);
    
    % plot
    subplot(2,2,1); hold on;
    plot(xDistance,nDistance(i,:),colors{i});
    subplot(2,2,2); hold on;
    plot(xDuration,nDuration(i,:),colors{i});
    subplot(2,2,3); hold on;
    plot(xIsi,nIsiBefore(i,:),colors{i});
    subplot(2,2,4); hold on;
    plot(xIsi,nIsiAfter(i,:),colors{i});
    
    % record
    subjects{i} = sprintf('subj %d',results{i}(1).subject);
end

% Get mean
nDistanceMean = mean(nDistance,1);
nDurationMean = mean(nDuration,1);
nIsiBeforeMean = mean(nIsiBefore,1);
nIsiAfterMean = mean(nIsiAfter,1);
subjects{nSubjects+1} = 'Mean';
subplot(2,2,1); hold on;
plot(xDistance,nDistanceMean,'k--','linewidth',2);
subplot(2,2,2); hold on;
plot(xDuration,nDurationMean,'k--','linewidth',2);
subplot(2,2,3); hold on;
plot(xIsi,nIsiBeforeMean,'k--','linewidth',2);
subplot(2,2,4); hold on;
plot(xIsi,nIsiAfterMean,'k--','linewidth',2);

% Label axes
subplot(2,2,1); 
xlabel('Saccade Distance (pixels)')
ylabel('# saccades')
subplot(2,2,2); 
xlabel('Saccade Duration (ms)')
ylabel('# saccades')
subplot(2,2,3); 
xlabel('isi before this saccade (ms)')
ylabel('# saccades')
subplot(2,2,4); 
xlabel('isi after this saccade (ms)')
ylabel('# saccades')    
legend(subjects)
MakeFigureTitle(sprintf('Saccade Class %s',num2str(class)));