function [avgW, avgTC, steW, steTC, multipliers] = MatchComponents(weights,timecourse)

% Adjust components and their timecourses so they all have the same sign.
% 
% [avgW, avgTC, steW, steTC, multipliers] = MatchComponents(weights,timecourse)
% 
% Created 3/13/14 by DJ based on GetGroupSvdResults.

nSubjects = numel(weights);
nComponents = size(weights{1},2);
nEventTypes = size(timecourse{1},1);
% Adjust components so all components have the "same sign"
isok = zeros(nSubjects,nComponents);
multipliers = ones(nSubjects,nComponents);
while ~all(isok(:))
    isok(:) = 1;
    avgW_temp = GetAverageWeights(weights);
    for i=1:nSubjects
        for j=1:nComponents
%             fprintf('S%s, c%d: %g\n',S(i).EEG.subject,j,weights{i}(iChans,j)'*avgW_temp(:,j));
            if weights{i}(:,j)'*avgW_temp(:,j)<0 % check whether individual and group are coherent
                fprintf('Reversing sign on subject %d, component %d...\n',i,j);
                multipliers(i,j) = -multipliers(i,j);
                weights{i}(:,j) = -weights{i}(:,j); % if not, reverse sign on individual's weights
                timecourse{i}(:,:,j) = -timecourse{i}(:,:,j); % and reverse sign on corresponding timecourse. 
                isok(i,j) = 0;
            end
        end
    end
end

% Get average component
[avgW,steW] = GetAverageWeights(weights);

% Get averages and standard errors
avgTC = mean(cat(4,timecourse{:}),4);
steTC = std(cat(4,timecourse{:}),[],4)/sqrt(nSubjects);

disp('Done!');



% --- HELPER FUNCTION --- %
% Get average and stderr of weights
function [avgW, steW] = GetAverageWeights(weights)

% Reformat to make taking mean easier
nSubjects = numel(weights);
nComponents = size(weights{1},2);
nElecs = size(weights{1},1);
wts = cell(1,nComponents);
for i=1:nSubjects
    for j=1:nComponents
        wts{j}(:,i) = weights{i}(:,j);
    end
end
% Calculate up avg and ste for each component
avgW = zeros(nElecs,nComponents);
steW = zeros(size(avgW));
for j=1:nComponents
    avgW(:,j) = mean(wts{j},2);
    steW(:,j) = std(wts{j},[],2)/sqrt(nSubjects);
end
