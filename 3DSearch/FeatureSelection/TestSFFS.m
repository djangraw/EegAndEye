function AUC = TestSFFS(class1,class2)

% Created 8/23/13 by DJ for one-time use.

truth = [zeros(1,size(class1,2)), ones(1,size(class2,2))];

y = [mean(class1,1), mean(class2, 1)];

AUC = rocarea(y,truth);