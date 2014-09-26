function Test3dsSffs_LOSubjO_script()

% Created 8/27/13 by DJ.


%% Get features of size X
nFeats = size(bigFeature{1},2);
isInSet = zeros(NumFeatComb,nFeats,numel(subjects));
for iSubj=1:numel(subjects)
    thisset = [];
    for i=1:NumFeatComb
%         if ~isempty(sffs.features{iSubj}{i});
            thisset = sffs.features{iSubj}{i};
%         end
        isInSet(i,thisset,iSubj) = 1;
    end
end


%% Get cross-validated AUC by leaving one subject out
bestfeats = cell(nFeats,numel(subjects));
AUC_xval = zeros(nFeats,numel(subjects));
firstmajority = zeros(nFeats,numel(subjects)); % to record any ties
for iSubj = 1:numel(subjects)   
    fprintf('====== Subject %d ======\n',subjects(iSubj));
    % Rank best features from training subjects
    trainsubj = [1:iSubj-1, iSubj+1:numel(subjects)];
    fracpicked = mean(isInSet(:,:,trainsubj),3);    
    for i=1:nFeats
        firstmajority(i,iSubj) = find(fracpicked(:,i)>0.5,1);
    end
    [~,featureorder] = sort(firstmajority(:,iSubj),'ascend');        
        
    % Record and 
    for nFeatWinners = 1:nFeats           
        bestfeats{nFeatWinners,iSubj} = featureorder(1:nFeatWinners);
        fprintf('Best %d features: [%s]\n',nFeatWinners,num2str(bestfeats{nFeatWinners,iSubj}'));
        % Use to test
        AUC_xval(nFeatWinners,iSubj) = ClassifierWrapper_forSFFS(bigFeature{iSubj}(:,bestfeats),bigFeature{iSubj}(:,bestfeats),EEGDATA{iSubj},YDATA{iSubj},offsets(iSubj));
    end        
end
