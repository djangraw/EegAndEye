%%% by Jun Wang EE@Columbia (jwang@ee.columbia.edu);
%%% Extended work of Graph Transduction via Alternating Minimization, ICML
%%% 2008

%%% Initial
%%% Oct. 21 2007 version by Jun Wang
%%% Propagate the highest score of F matrix;
%%% updating by pick up the most possible samples; This version is robust
%%% to weight matrix;
%%% Oct. 22 2008 test noisy and uneven cases
%%% Oct. 24 2008 revise to be more robust to noise.
%%% Revise the node term;
%%% May. 12, 2008
%%% Update July 12, 2009
%%% Updated July 16, 2013 by DJ - added labeled_ind output


function [F,labeled_ind]=GTAM_eegwronglabels(labeled_ind,A,W,IS,options)


%%% Rerank EEG test results
%%% For the experiemnts on ACM MM 2009 Paper "Brain State Decoding for Rapid Image Retrieval"
%%%                           Here using a standard way. You can construct
%%%                           your own graph
%%% Input:
%%%       --- labeled_ind      Index of EEG labels
%%%       --- A                graph gradient matrix (Equation 9 in ICML 08 paper) = graph.gradient
%%%       --- W                weight matrix = graph.weight
%%%       --- IS               propagation matrix = graph.propm     
%%%       --- options          how many iterations for label self tuning?; (only field rmvnum is used)

%%% Output:
%%%       --- F                re-ranking results;
%%%       --- labeled_ind      EEG labels after self-tuning
%%%
%%% by J. Wang (jwang@ee.columbia.edu)
%%% Last update July 16, 2013

%%%only for debugging;
%load imgname_list.mat;
data_num=length(A); % number of images
class_num=1;

%%% coefficient for local label fitting

%%% iteration counter
iter_num=0;
%%% predicted labels;
predict_label=zeros(1,data_num);
predict_label(labeled_ind)=1;
%%% unlabeled index
unlabeled_ind=setdiff([1:data_num],labeled_ind);

%%% Initial given labels
Y=zeros(data_num,class_num);
for i=1:length(labeled_ind)
    ii=labeled_ind(i);
    Y(ii)=1;     
end
originalY=Y;

%%% Initial node regularizer
subW=W(labeled_ind,labeled_ind);
subD=sum(subW);
V=zeros(data_num,data_num);
for i=1:length(labeled_ind)
    ii=labeled_ind(i);
    %V(ii,ii)=d(ii)/sum(d(predict_label==1));     
    V(ii,ii)=max([0.001 subD(i)/sum(subD)]);
end
    

F=zeros(data_num,1);
while (isempty(unlabeled_ind)~=1)&&(iter_num<=options.rmvnum)   
%%%   Partial differential to VY    
    normalizedY=V*Y;
    DeltaQ=A*normalizedY;
    
    %%% Remove low-confidence label;
    DeltaY=-DeltaQ(labeled_ind,:);
    [a, b]=max(DeltaY);
    predict_label(labeled_ind(b))=0;
    %%% for debugging;
    %fprintf('remove %s \n',imgname_list{labeled_ind(b)});
    labeled_ind(b)=[];
    unlabeled_ind=setdiff([1:data_num],labeled_ind);    
    
    %%%refine the label weights
    subW=W(labeled_ind,labeled_ind);
    subD=sum(subW);
    clear V;
    V=zeros(data_num,data_num);    
    for i=1:length(labeled_ind)
        ii=labeled_ind(i);
        %V(ii,ii)=d(ii)/sum(d(predict_label==predict_label(ii)));   
        V(ii,ii)=max([0.001 subD(i)/sum(subD)]);
    end    
    iter_num=iter_num+1;
    Y=zeros(data_num,class_num);
    for i=1:length(labeled_ind)
        ii=labeled_ind(i);
        Y(ii)=1;     
    end

    %%% Add high-confidence label;
    normalizedY=V*Y;
    DeltaQ=A*normalizedY;    
    DeltaY=-DeltaQ(unlabeled_ind,:);
    [a, b]=max(DeltaY);
    labeled_ind=[labeled_ind,unlabeled_ind(b)];
    predict_label(unlabeled_ind(b))=1;    
    %%% for debugging;
    %fprintf('add %s \n', imgname_list{unlabeled_ind(b)});   
    unlabeled_ind=setdiff([1:data_num],labeled_ind);
  
   
    subW=W(labeled_ind,labeled_ind);
    subD=sum(subW);    
    %%%refine the label weights
    clear V;
    V=zeros(data_num,data_num);    
    for i=1:length(labeled_ind)
        ii=labeled_ind(i);
        %V(ii,ii)=d(ii)/sum(d(predict_label==predict_label(ii)));  
        V(ii,ii)=max([0.001 subD(i)/sum(subD)]);
    end    
    Y=zeros(data_num,class_num);
    for i=1:length(labeled_ind)
        ii=labeled_ind(i);
        Y(ii)=1;     
    end    
end

% while (isempty(unlabeled_ind)~=1)&&(iter_num<=5)   
%     %%% Add high-confidence label;
%     normalizedY=V*Y;
%     DeltaQ=A*normalizedY;    
%     DeltaY=-DeltaQ(unlabeled_ind,:);
%     [a, b]=max(DeltaY);
%     labeled_ind=[labeled_ind,unlabeled_ind(b)];
%     predict_label(unlabeled_ind(b))=1;    
%     fprintf('add %s \n', imgname_list{unlabeled_ind(b)});   
%     unlabeled_ind=setdiff([1:data_num],labeled_ind);
%   
%    
%     subW=W(labeled_ind,labeled_ind);
%     subD=sum(subW);    
%     %%%refine the label weights
%     clear V;
%     V=zeros(data_num,data_num);    
%     for i=1:length(labeled_ind)
%         ii=labeled_ind(i);
%         %V(ii,ii)=d(ii)/sum(d(predict_label==predict_label(ii)));  
%         V(ii,ii)=subD(i)/sum(subD);
%     end    
%     Y=zeros(data_num,class_num);
%     for i=1:length(labeled_ind)
%         ii=labeled_ind(i);
%         Y(ii)=1;     
%     end    
% end

Y=zeros(data_num,class_num);
for i=1:length(labeled_ind)
    ii=labeled_ind(i);
    Y(ii)=1;         
end

clear V;
V=zeros(data_num,data_num);

subW=W(labeled_ind,labeled_ind);
subD=sum(subW);
for i=1:length(labeled_ind)
    ii=labeled_ind(i);
    %V(ii,ii)=d(ii)/sum(d(predict_label==predict_label(ii)));  
    V(ii,ii)=max([0.001 subD(i)/sum(subD)]);
end   
clear W;
F=IS*V*Y;
%F=F.*Y;



