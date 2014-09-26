function [cLBest,maxJ]=SequentialForwardSelection(class1,class2,CostFunction,NumFeatComb,otherInputs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
%   [cLBest,maxJ]=SequentialForwardSelection(class1,class2,CostFunction,NumFeatComb,otherInputs);
%   Feature vector selection by means of the Sequential Forward Selection
%   technique, given the desired number of features in the best combination.
%
% INPUT ARGUMENTS:
%   class1:         matrix of data for the first class, one pattern per column.
%   class2:         matrix of data for the second class, one pattern per column.
%   CostFunction:   class separability measure.
%   NumFeatComb:    desired number of features in best combination.
%   otherInputs:    cell array of other arguments used as input to cost fn.
%
% OUTPUT ARGUMENTS:
%   cLBest:         selected feature subset. Vector of row indices.
%   maxJ:           class separabilty measure.
%
% (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
%
% Updated 8/23/13 by DJ - added otherInputs input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NofFeatures=size(class1,1);
cLBest=[];
k=1;
while k<=NumFeatComb
    maxJ=0;
    for i=1:NofFeatures
        if ~any(cLBest==i) % DJ change
            combi=[cLBest i];
        else continue;
        end
        eval(['[J]=' CostFunction '(class1(combi,:),class2(combi,:),otherInputs{:});']);
        if J>maxJ
            maxJ=J;
            sofar=combi;
        end
    end
    cLBest=sort(sofar,'ascend');
    k=k+1;
end
