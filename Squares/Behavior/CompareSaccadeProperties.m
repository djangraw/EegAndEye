function CompareSaccadeProperties(results,classes)

% Plot histograms of saccade distance, duration, and time to previous or
% next saccade for multiple classes, averaging across subjects
%
% CompareSaccadeProperties(results,classes)
%
% INPUTS:
% -results is a cell array in which each cell contains a squares dataset 
% imported by import_squares_data.m for one subject.
% -classes is a vector of event types (see UseEventRule.m for details) 
% indicating the saccade type you want to analyze.
%
% Created 7/16/12 by DJ based on GetSaccadeProperties
% Updated 7/23/12 by DJ - added SquareDuration measure

% Declare constants
nSubjects = length(results);
nClasses = length(classes);
colors = {'r' 'g' 'b' 'c' 'm' 'y' 'r--' 'g--' 'b--' 'c--' 'm--' 'y--'};
% Declare histogram bins
xDistance = 5:10:300;
xDuration = 2:4:100;
xIsi = 10:20:1000;

% Set up loop
nDistance = zeros(nClasses,length(xDistance),nSubjects);
nDuration = zeros(nClasses,length(xDuration),nSubjects);
nIsiBefore = zeros(nClasses,length(xIsi),nSubjects);
nIsiAfter = zeros(nClasses,length(xIsi),nSubjects);
nEvents = zeros(nClasses,nSubjects);
clf; % clear current figure

% Main loop
for i=1:nSubjects
    % Find relevant saccades
%     nSessions = length(results{i});        
%     s = combine_structs([results{i}.saccade],[],'columns');
    
    for j=1:nClasses    
        % Get event indices
        [~,~,~,iEvents] = UseEventRule(results{i},classes{j});
        % Get event info
        nSessions = numel(results{i});
        [dist,dur,squaredurA,fixdurA] = deal(cell(1,nSessions));
        for k=1:nSessions
            s = results{i}(k).saccade;
            dist{k} = sqrt((s.end_position(iEvents{k},1)-s.start_position(iEvents{k},1)).^2 + ...
            (s.end_position(iEvents{k},2)-s.start_position(iEvents{k},2)).^2);
            dur{k} = s.end_time(iEvents{k}) - s.start_time(iEvents{k});
            % Get duration of fixating on this square
            squaredur = nan(size(s.end_time));
            iChangeSquare = find(diff(s.squarenum)~=0)+1; % saccades when subject moved away from current square
            squaredur(iChangeSquare(1:end-1)) = s.start_time(iChangeSquare(2:end)) - s.end_time(iChangeSquare(1:end-1));
            squaredurA{k} = squaredur(iEvents{k});
            % Get fixation duration
            fixdur = s.start_time(2:end) - s.end_time(1:end-1);
%             fixdurB{k} = fixdur(iEvents{k}-1); % duration of fixation before this saccade
            fixdurA{k} = fixdur(iEvents{k}); % duration of fixation after this saccade
        end
        % Concatenate across sessions
        distance = cat(1,dist{:});
        duration = cat(1,dur{:});
        squaredur_after = cat(1,squaredurA{:});
%         fixdur_before = cat(1,fixdurB{:});
        fixdur_after = cat(1,fixdurA{:});
        nEvents(j,i) = length(distance);
        % Run histograms
        nDistance(j,:,i) = hist(distance,xDistance)/nEvents(j,i);
        nDuration(j,:,i) = hist(duration,xDuration)/nEvents(j,i);
        nSDAfter(j,:,i) = hist(squaredur_after,xIsi)/nEvents(j,i);
%         nFDBefore(j,:,i) = hist(fixdur_before,xIsi)/nEvents(j,i);
        nFDAfter(j,:,i) = hist(fixdur_after,xIsi)/nEvents(j,i);
    end
      
end

% Get mean
nDistanceMean = mean(nDistance,3);
nDurationMean = mean(nDuration,3);
nSDAfterMean = mean(nSDAfter,3);
% nFDBeforeMean = mean(nFDBefore,3);
nFDAfterMean = mean(nFDAfter,3);
nEventsMean = mean(nEvents,2);
% Get stddev
nDistanceStd = std(nDistance,0,3)/sqrt(nSubjects);
nDurationStd = std(nDuration,0,3)/sqrt(nSubjects);
nSDAfterStd = std(nSDAfter,0,3)/sqrt(nSubjects);
% nFDBeforeStd = std(nFDBefore,0,3)/sqrt(nSubjects);
nFDAfterStd = std(nFDAfter,0,3)/sqrt(nSubjects);



for j=1:nClasses
    subplot(2,2,1); hold on;
%     JackKnife(xDistance,nDistanceMean(j,:)*100,nDistanceStd(j,:)*100,colors{j},colors{j});
    plot(xDistance,nDistanceMean(j,:)*100,colors{j},'linewidth',2);
    subplot(2,2,2); hold on;
%     JackKnife(xDuration,nDurationMean(j,:)*100,nDurationStd(j,:)*100,colors{j},colors{j});
    plot(xDuration,nDurationMean(j,:)*100,colors{j},'linewidth',2);
    subplot(2,2,3); hold on;
%     JackKnife(xIsi,nSDAfterMean(j,:)*100,nSDAfterStd(j,:)*100,colors{j},colors{j});
%     plot(xIsi,nFDBeforeMean(j,:)*100,colors{j},'linewidth',2);
    plot(xIsi,nSDAfterMean(j,:)*100,colors{j},'linewidth',2);
    subplot(2,2,4); hold on;
%     JackKnife(xIsi,nFDAfterMean(j,:)*100,nFDAfterStd(j,:)*100,colors{j},colors{j});
    plot(xIsi,nFDAfterMean(j,:)*100,colors{j},'linewidth',2);
end
    
% Label axes
subplot(2,2,1); 
ylim([0,max(get(gca,'ylim'))]);
xlabel('Saccade Distance (pixels)')
ylabel('% of saccades')
subplot(2,2,2); 
ylim([0,max(get(gca,'ylim'))]);
xlabel('Saccade Duration (ms)')
ylabel('% of saccades')
subplot(2,2,3); 
ylim([0,max(get(gca,'ylim'))]);
xlabel('duration of fixation on this square (ms)')
% xlabel('duration of fixation before this saccade (ms)')
ylabel('% of saccades')
subplot(2,2,4); 
ylim([0,max(get(gca,'ylim'))]);
xlabel('duration of fixation after this saccade (ms)')
ylabel('% of saccades')
% Make legend
classnames = cell(1,nClasses);
for j=1:nClasses
    classnames{j} = sprintf('%s (n=%d)',classes{j},round(nEventsMean(j)));
end
MakeLegend(colors(1:nClasses),classnames,repmat(2,1,nClasses));
MakeFigureTitle(sprintf('Saccade Properties across %d subjects',nSubjects));