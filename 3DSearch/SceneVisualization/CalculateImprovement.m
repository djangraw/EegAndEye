function [stats,TagTour,plotInfo] = CalculateImprovement(subject,sessions,levelname,eegPredictedTargets,usegridconstraints,nSensitivity)

% Calculates targets seen and distances traveled.
%
% [stats,TagTour,plotInfo] = CalculateImprovement(subject,sessions,levelname,eegPredictedTargets,usegridconstraints,nSensitivity)
%
% INPUTS:
% -subject is the subject number and sessions is a vector of sessions
% numbers of 3DSearch data files as imported with Import_3DS_Data_v3.
% -levelname is a string that indicates an image filename in the current
% path.
% - eegPredictedTargets is a vector of the object indices (in the combined
% dataset) predicted to be targets based on the EEG.
% -usegridconstraints is a boolean indicating whether you want the
% traveling salesman algorithm to snap all the object and routes to the
% grid (otherwise it can travel diagonally). [default: false]
% -nSensitivity is a scalar indicating how many times you want to perform
% "self-tuning" (throwing out one predicted target that doesn't match and
% adding another that does). [Default: 0]. fracSensitivity may also be used
% (see RerankObjectsWithTag for details).
%
% OUTPUT:
% -stats is a struct with some information from the analysis.
% -TagTour is an Nx2 matrix, where N is the number of stops on the
% traveling salesman route suggested by this program.
% -plotInfo is a struct with info that can be used to plot the results.
%
% Created 3/9/12 by DJ based on ImageAllSessions_Bionav.
% Updated 8/20/12 by DJ - added stats output.
% Updated 12/13/12 by DJ - added TagTour output.
% Updated 5/13/13 by DJ - added optional nSensitivity input.
% Updated 6/7/13 by DJ - added plotInfo output, fixed nObjSeenPerSession.
% Updated 7/12/13 by DJ - switched to guess of 95% point as cutoff.
% Updated 7/17/13 by DJ - switched to MoG fit to find cutoff.

if nargin<5 || isempty(usegridconstraints)
    usegridconstraints = false;
end
if nargin<6 || isempty(nSensitivity)
    nSensitivity = 0;
end
% Set up
nCategories = 4;

% load camera path
levelinfo = load([levelname(1:end-4) '.mat']); % load variables campoints and sessionOffset

if strcmp(levelname(1:end-4),'GridHuge')
    nPointsPerSession = 23;
    nObjSeenPerSession = 20;
else
    nPointsPerSession = 21;
end

% load object info
[objects, objnames, objlocs, objtimes] = GetObjectList(subject,sessions);

% main loop
camdist = 0;
cameraPath = cell(1,numel(sessions));
for i=1:numel(sessions)
    % initialize
    load(sprintf('3DS-%d-%d',subject,sessions(i)));
    nObjPerSession = length(x.objects);
   
    % get spatial offset for this session
    thisSessionOffset = levelinfo.sessionOffset * (i-1); % this offset specific to level SnakeHuge

    % add distance traveled to record
    nObjectsSeen = cumsum(levelinfo.isObjPoint((x.start_point+2):end));
    iLastPoint = find(nObjectsSeen==nObjSeenPerSession);
    campath = levelinfo.campoints(x.start_point+1+(1:iLastPoint),:); 
    pathlength = sum(sqrt(sum(diff(campath,1).^2,2)));        
    camdist = camdist+pathlength;
    % save camera path
    cameraPath{i} = campath + repmat(thisSessionOffset,length(campath),1);
    
    % shift objects
    iThisSession = (i-1)*nObjPerSession + (1:nObjPerSession); % indices in obj___ for this session's objects
    objlocs(iThisSession,:) = objlocs(iThisSession,:) + repmat(thisSessionOffset,nObjPerSession,1);    

end

% Snap object locations to grid
if usegridconstraints
    objlocs(:,1) = round(objlocs(:,1)/15)*15;
    objlocs(:,2) = round(objlocs(:,2)/20)*20;
end


% categorize points according to true label
targets = strcmp('TargetObject',{objects(:).tag});
distractors = strcmp('DistractorObject',{objects(:).tag});

% further separate points that were seen
wasViewed = ~isnan(objtimes);

% Show ideal result of EEG classifier
if nargin<4
    eegPredictedTargets = find(targets & wasViewed); % default: perfect classification
end

% Get TAG predicted-target objects
[newRanking,outScore,isSelfTunedTarget] = RerankObjectsWithTag(objnames, eegPredictedTargets,nSensitivity);
% % Find top 1/nCategories of objects
% rankCutoff = round(length(newRanking)/nCategories);
% iTagTargets = newRanking(1:rankCutoff);

% % Find "excessively large" TAG scores
% % TEMP_PlotTargetFracPrediction;
% fivepct = outScore(round(numel(outScore)*0.95));
% fiftypct = median(outScore);
% ninetyfivepct = 2*fiftypct-fivepct;
% iTagTargets = newRanking(outScore>ninetyfivepct);

%% Fit to MoG (Gaussian mixture model)
obj = gmdistribution.fit(outScore',2);
% Find intersection point
% compute solutions when the sigmas are different
A1 = obj.PComponents(1);
A2 = obj.PComponents(2);
mu1 = obj.mu(1);
mu2 = obj.mu(2);
sigma1 = sqrt(obj.Sigma(1));
sigma2 = sqrt(obj.Sigma(2));
SIGMA1 = sigma1 .* sigma1;
SIGMA2 = sigma2 .* sigma2;

aux1 = mu2 .* SIGMA1 - mu1 .* SIGMA2;
aux2 = sigma1 .* sigma2 .* sqrt((...
    (mu1 - mu2).^2 + 2 * (SIGMA2 - SIGMA1) .* log(A1.*sigma2 ./ (A2.*sigma1)) ...
    ));
aux3 = 1./(SIGMA1 - SIGMA2);

% compute intersection points
point(1) = (aux1 - aux2) .* aux3;
point(2) = (aux1 + aux2) .* aux3;
fprintf('TAG cutoff is %f\n',max(point))
iTagTargets = newRanking(outScore>max(point));
%%

% Get binary bits (useful below)
isEegTarget = zeros(size(targets));
isEegTarget(eegPredictedTargets) = 1;
isTagTarget = zeros(size(targets));
isTagTarget(iTagTargets) = 1;
iSelfTunedTargets = newRanking(isSelfTunedTarget==1);
isSelfTunedTarget = zeros(size(targets));
isSelfTunedTarget(iSelfTunedTargets) = 1;

% Get traveling salesman route
TagLocations = objlocs(iTagTargets,:);
[TagTour,~,TagTourDist] = solveTSP(TagLocations,0,usegridconstraints);
[AllTour,~,AllTourDist] = solveTSP(objlocs,0,usegridconstraints);

% Display results
fprintf('--------------------\n')
fprintf('Percent correct -- EEG: %0.1f, TAG: %0.1f\n',...
    sum(isEegTarget & targets)/sum(isEegTarget)*100,...
    sum(isTagTarget & targets)/sum(isTagTarget)*100);
fprintf('Percent of targets found -- EEG: %0.1f, TAG: %0.1f\n',...
    sum(isEegTarget & targets)/sum(targets)*100,...
    sum(isTagTarget & targets)/sum(targets)*100);
fprintf('%% of TAG pts that are true targets: %0.1f\n', ...
    sum(isTagTarget & targets)/sum(isTagTarget)*100);
fprintf('%% of true targets that are TAG pts: %0.1f\n', ...
    sum(isTagTarget & targets)/sum(targets)*100);
fprintf('--------------------\n')
fprintf('Total objects in scene: %d\n', ...
    numel(objects));
fprintf('Targets seen during training: %d/%d\n', ...
    sum(targets & wasViewed), sum(wasViewed));
fprintf('Targets seen during Tag traveling salesman: %d/%d\n', ...
    sum(targets & isTagTarget), sum(isTagTarget));
fprintf('New targets seen during Tag traveling salesman: %d\n', ...
    sum(targets & isTagTarget & ~wasViewed));
fprintf('--------------------\n')
fprintf('Distance traveled during training: %0.1fm\n',camdist);
fprintf('All objects traveling salesman distance: %0.1fm\n', AllTourDist);
fprintf('Tag traveling salesman distance: %0.1fm\n', TagTourDist);   
fprintf('Savings: travel %0.1f%% of the distance to see %0.1f%% of the targets.\n',...
    TagTourDist/AllTourDist*100, sum(isTagTarget & targets)/sum(targets)*100);
fprintf('--------------------\n')    

stats.precision(1) = sum(isEegTarget & targets)/sum(isEegTarget)*100;
% stats.precision(2) = sum(isSelfTunedTarget & targets)/sum(isSelfTunedTarget)*100;
stats.precision(2) = sum(isTagTarget & targets)/sum(isTagTarget)*100;
stats.pctFound(1) = sum(isEegTarget & targets)/sum(targets)*100;
% stats.pctFound(2) = sum(isSelfTunedTarget & targets)/sum(targets)*100;
stats.pctFound(2) = sum(isTagTarget & targets)/sum(targets)*100;
stats.pctDistance = TagTourDist/AllTourDist*100;

% Make plotInfo struct for plotting
plotInfo.objlocs = objlocs;
plotInfo.iTargets = targets;
plotInfo.iEegTargets = eegPredictedTargets;
plotInfo.iSelfTunedTargets = iSelfTunedTargets;
plotInfo.iTagTargets = iTagTargets;
% tagScore(newRanking) = length(newRanking):-1:1; % first-ranked obj has highest score
% plotInfo.tagScore = tagScore;
plotInfo.tagScore(newRanking) = outScore;
plotInfo.tagRank(newRanking) = 1:length(newRanking);
plotInfo.cameraPath = cameraPath;