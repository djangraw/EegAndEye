function ImageAllSessions(subject,sessions,levelname,offset)

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

% display level
ImageLevel(levelname,offset);

% load camera path
load([levelname(1:end-4) '.mat']) % load variable campoints

% main loop
for i=1:numel(sessions)
    % initialize
    target_points = [];
    distractor_points = [];
    load(sprintf('3DS-%d-%d.mat',subject,sessions(i))); % get x struct
    % add each object to target or distractor point list
    for j=1:numel(x.objects)
        if strcmp(x.objects(j).tag, 'DistractorObject')
            distractor_points = [distractor_points; ...
                x.objects(j).createposition(1), x.objects(j).createposition(3)];
        elseif strcmp(x.objects(j).tag,'TargetObject')
            target_points = [target_points; ...
                x.objects(j).createposition(1), x.objects(j).createposition(3)];
        end
    end
    % get camera path
    campath = campoints(x.start_point+1+(1:20),:);
    
    % plot with offset
    thisSessionOffset = [80 100] * (i-1); % this offset specific to level SnakeHuge
    ImagePoints(target_points,thisSessionOffset,'r','.'); % red targets
    ImagePoints(distractor_points,thisSessionOffset,'g','.'); % green distractors
    ImagePoints(campath,thisSessionOffset,'c','c--'); % blue camera path
end

