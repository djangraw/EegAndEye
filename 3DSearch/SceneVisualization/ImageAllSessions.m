function ImageAllSessions(subject,sessions,levelname,offset,eegPredictedTargets)

% Displays an image of the indicated level with the camera position, target
% and distractor objects superimposed.
%
% ImageAllSessions(subject,sessions,levelname,offset)
%
% INPUTS:
% -subject is the subject number and sessions is a vector of sessions
% numbers of 3DSearch data files as imported with Import_3DS_Data_v3.
% -levelname is a string that indicates an image filename in the current
% path.
% - offset is the x,y position of the origin in the image.  It will be
% subtracted from the x and y axes when plotting the image. [default: 0 0]
%
% Created 4/18/11 by DJ.

% handle inputs
if nargin<4
    offset = [0 0];
end

nCategories = 4;

% display level
ImageLevel(levelname,offset);

% load camera path
levelinfo = load([levelname(1:end-4) '.mat']); % fields campoints and sessionOffset

% load object info
[objects, objnames, objlocs, objtimes] = GetObjectList(subject,sessions);

% main loop
for i=1:numel(sessions)
    % initialize
    load(sprintf('3DS-%d-%d',subject,sessions(i)));
    nObjPerSession = length(x.objects);
   
    % get spatial offset for this session
    thisSessionOffset = levelinfo.sessionOffset * (i-1); % this offset specific to level SnakeHuge

    % plot camera path with offset
    campath = levelinfo.campoints(x.start_point+1+(1:21),:); 
    ImagePoints(campath,thisSessionOffset,'c','c--'); % blue camera path    
    
    % shift objects
    iThisSession = (i-1)*nObjPerSession + (1:nObjPerSession); % indices in obj___ for this session's objects
    objlocs(iThisSession,:) = objlocs(iThisSession,:) + repmat(thisSessionOffset,nObjPerSession,1);    

end

% categorize points according to true label
targets = strcmp('TargetObject',{objects(:).tag});
distractors = strcmp('DistractorObject',{objects(:).tag});

% further separate points that were seen
wasViewed = ~isnan(objtimes);

% plot true labeled points
ImagePoints(objlocs(targets,:),[0 0],'r','.'); % red targets
ImagePoints(objlocs(distractors,:),[0 0],'b','.'); % green distractors

% Show ideal result of EEG classifier
if nargin<5
    eegPredictedTargets = find(targets & wasViewed); % default: perfect classification
end
ImagePoints(objlocs(eegPredictedTargets,:),[0 0],'y','o',10); % yellow targets

% Show TAG's leveraging of ideal classifier
newRanking = RerankObjectsWithTag(objnames, eegPredictedTargets);
iTagTargets = newRanking(1:end/nCategories);
ImagePoints(objlocs(iTagTargets,:),[0 0],'g','o',15); % green targets

% Annotate plot
title(sprintf('Subject %d, level %s', subject, levelname(1:end-4)));
MakeLegend({'c--','r.','b.','yo','go'},...
    {'subject route', 'targets','distractors','EEG predicted targets', 'TAG predicted targets'})
