function [outputOrder,outputScore,isSelfTunedTarget] = RerankObjectsWithTag(objectList,iTargets,nSensitivity)

% Use a TAG graph to re-rank the given objects from a 3DSearch task.
%
% outputOrder = RerankObjectsWithTag(objectList,iTargets,nSensitivity)
% outputOrder = RerankObjectsWithTag(objectList,iTargets,fracSensitivity)
%
% INPUTS:
% -objectList is a 1xn array of strings, each of format
% '<category>-<number>'.  This indicates the category and number of the
% CalTech-101 image used to create this prefab.
% -iTargets is an array of objectIndices that were indicated as targets by
% the EEG classifier.
% -nSensitivity is a scalar indicating how many times you want to perform
% "self-tuning" (throwing out one predicted target that doesn't match and
% adding another that does). [Default: 0].
% -fracSensitivity is a scalar <1 indicating what fraction of the targets
% should be used as the number of times to perform "self-tuning". (If you
% want to perform 100% self-tuning, try 0.9999999).
%
% OUTPUTS:
% -outputOrder is the integers 1:n in the order of their tag rankings.
% -outputScore is the score output by TAG for the corresponding images.
% -iSelfTunedTargets(i) indicates whether object i was a predicted target
% after self-tuning.
%
% Created 4/18/11 by DJ.
% Updated 1/15/13 by DJ - added useLite option (only includes 4 categories
% used in 3DSearch <v7p1).
% Updated 5/13/13 by DJ - added optional nSensitivity input.
% Updated 5/16/13 by DJ - added fracSensitivity option.
% Updated 7/12/13 by DJ - added outputScore output.
% Updated 7/16/13 by DJ - added iSelfTunedTargets output.


%% Set up
% indicate amount of self-tuning
if nargin<3 || isempty(nSensitivity)
    nSensitivity = 0;
end

% interpret fracSensitivity option
if nSensitivity < 1
    fracSensitivity = nSensitivity;
    nSensitivity = round(numel(iTargets)*fracSensitivity);    
end
fprintf('--- TAG: Self-tuning will be applied %d times.\n',nSensitivity);

%%% load image list & TAG graph
graphtype = 'tiny'; % tiny: only include images in Unity. lite: only include images from 4 categories.
switch graphtype
    case 'tiny'
        fprintf('Loading FileList_tiny and graph_3_tiny...\n');
        load FileList_tiny
        load graph_3_tiny
    case 'lite'
        fprintf('Loading FileList_lite and graph_3_lite...\n');
        load FileList_lite
        load graph_3_lite
    otherwise
        graphnum = 3;
        fprintf('Loading FileList and graph_%d...\n', graphnum);
        load FileList    
        load(['graph_' num2str(graphnum) '.mat'])
end
%% Get indices in TAG graph
disp('Getting TAG indices of objects...')
objectInd = zeros(1,numel(objectList));
targetInd = [];
for i=1:numel(objectList)
    values = textscan(objectList{i},'%s %d','Delimiter','-');
    catname = values{1}{1};
    imnum = values{2};
    objectInd(i) = strmatch(sprintf('%s_image_%04.f.jpg',catname,imnum),FileList);
    if ismember(i,iTargets) % add to target list if appropriate
        targetInd = [targetInd objectInd(i)];
    end    
end

%% Run TAG
disp('Re-ranking all images with TAG...')
%%% Set options
options.rmvnum=nSensitivity;     %%% Self tuning for 10 iterations;

%%% it takes a while to run reranking
[new_score, iSelfTunedTargets] = GTAM_eegwronglabels(targetInd,graph.gradient,graph.weight,graph.propm,options);

isSST_tagout = zeros(size(new_score));
isSST_tagout(iSelfTunedTargets) = 1;

%% Re-rank input images accordingly
disp('Re-ranking objects...')
[re_rank_score, re_rank_ind]=sort(new_score,'descend');
iReranked = re_rank_ind(ismember(re_rank_ind,objectInd));
newscore = re_rank_score(ismember(re_rank_ind,objectInd));

isSST_rerankscore = isSST_tagout(re_rank_ind);
isSST_newscore = isSST_rerankscore(ismember(re_rank_ind,objectInd));
outputOrder = [];
outputScore = [];
isSelfTunedTarget = [];
for i=1:numel(iReranked)
    outputOrder = [outputOrder find(objectInd==iReranked(i))];
    outputScore = [outputScore repmat(newscore(i),1,sum(objectInd==iReranked(i)))];
    isSelfTunedTarget = [isSelfTunedTarget repmat(isSST_newscore(i),1,sum(objectInd==iReranked(i)))];
end

disp('Success!')