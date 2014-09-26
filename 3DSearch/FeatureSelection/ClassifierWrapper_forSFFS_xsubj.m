function [AUC,v] = ClassifierWrapper_forSFFS_xsubj(bigFeature,~,truth)

% Created 8/27/13 by DJ.
% Updated 8/28/13 by DJ - added v output

% Make bigFeature into one cell per subject
nSubjects = size(bigFeature,2);
features = cell(1,nSubjects);
for iSubj=1:nSubjects
    features{iSubj} = cat(2,bigFeature{:,iSubj});
end

% classify
[subjAUC,v] = ClassifyWithOcularFeatures(truth,features);
    
AUC = mean(subjAUC);