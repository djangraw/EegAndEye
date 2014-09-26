function [cLbest,maxJ]=SequentialForwardFloatingSelection(class1,class2,CostFunction,NumFeatComb,otherInputs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
%  [cLbest,maxJ]=SequentialForwardFloatingSelection(class1,class2,CostFunction,NumFeatComb);
%  Feature vector selection by means of the Sequential Forward Floating
%  Selection technique, given the desired number of features in the best combination.
%
% INPUT ARGUMENTS:
%   class1:         matrix of data for the first class, one pattern per column.
%   class2:         matrix of data for the second class, one pattern per column.
%   CostFunction:   class separability measure.
%   NumFeatComb:    desired number of features in best combination.
%   otherInputs:    cell array of other arguments used as input to cost fn.
%
% OUTPUT ARGUMENTS:
%   cLbest:         selected feature subset. Vector of row indices.
%   maxJ:           value of the class separabilty measure.
%
% (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
%
% Updated 8/23/13 by DJ - added otherInputs input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off');
m=size(class1,1);   % total # features
k=2;                % current # features
% Initialization
[X{k},C{k}]=SequentialForwardSelection(class1,class2,CostFunction,2,otherInputs); % best 2 features and classification score


while k<=NumFeatComb
    % Print current feature set
    fprintf('Set = [%s], cost = %g\n',num2str(X{k}),C{k})
    
    % Step I:Add    
    Ct=[]; % classification score of each option
    Y{m-k}=setdiff(1:m,X{k}); % features NOT in current set
    for i=1:length(Y{m-k})
        t=[X{k} Y{m-k}(i)]; % add each out feature to set
        Ct=[Ct eval([CostFunction '(class1(t,:),class2(t,:),otherInputs{:});'])]; % get classification score
    end
    [J_of_x,ind]=max(Ct); % Find best one
    the_x=Y{m-k}(ind); % Find additional feature
    X{k+1}=[X{k} the_x]; % Add to feature set
    
    % Step II:Test
    Ct=[];
    for i=1:length(X{k+1})
        t=setdiff(X{k+1}, X{k+1}(i)); % Remove each feature from set
        Ct=[Ct eval([CostFunction '(class1(t,:),class2(t,:),otherInputs{:});'])]; % get classification score
    end
    [J,r]=max(Ct); % Find max cost and feature index
    xr=X{k+1}(r); % Find feature to be removed
    if r==k+1 % If it's the one you just added
        C{k+1}=eval([CostFunction '(class1(X{k+1},:),class2(X{k+1},:),otherInputs{:});']); % Get classification score from adding it
        k=k+1; % Approve the addition
        continue;
    end
    if r~=k+1 &  J< C{k} % if removing something else wouldn't help classification score
        continue; % Stay with k, don't add new feature
    end
    if k==2
        X{k}=setdiff(X{k+1},xr);
        C{k}=J;
        continue;
    end
    
    % Step III: Exclusion
    flag=1;
    while flag
%         X_hat{k}=setdiff(X{k+1},xr); % MATLAB forum suggested removing this 
        X_hat{k} = X{k}; % DJ added this... is it right?
        Ct=[];
        for i=1:length(X_hat{k})
            t=setdiff(X_hat{k}, X_hat{k}(i));
            Ct=[Ct eval([CostFunction '(class1(t,:),class2(t,:),otherInputs{:});'])];
        end
        [J,s]=max(Ct);
        xs=X_hat{k}(s);
        if J<C{k-1}
            X{k}=X_hat{k};
            C{k}=eval([CostFunction '(class1(X{k},:),class2(X{k},:),otherInputs{:});']);
            flag=0;
            fprintf('Set = [%s], cost = %g\n',num2str(X{k}),C{k})
            break;
        end
        X_hat{k-1}=setdiff(X_hat{k},xs);
        k=k-1;
        if k==2
            X{k}=X_hat{k};
            C{k}=eval([CostFunction '(class1(X_hat{k},:),class2(X_hat{k},:),otherInputs{:});']);
            flag=0;
        end
        fprintf('Set = [%s], cost = %g\n',num2str(X{k}),C{k})
    end
    if flag==0
        continue;
    end
end

if k>NumFeatComb
    k=k-1;
end
cLbest=sort(X{k},'ascend');
maxJ=C{k};

