%=========================================================================%
%%% Iterative shrinkage algorithm for l1-regularized logistic regresssion
%%% Authors: Jianing Shi, Wotao Yin, Paul Sajda and Stanley Osher
%%% Contact: Jianing Shi   js2615@columbia.edu

%%% This algorithm solves l1-regularized logistic regression 
%%% Using an iterative shrinkage algorithm
%%% Given data X and b
%%% Seek the following optimization
%%% {w,v} = arg min L_log(w,v) + lambda ||w||_L1
%%% X : m*n
%%% b : m*1
%%% w : n*1
%%% v : 1*1
%%% lambda : 1*1

%%% (Copyright)  Jianing Shi, Wotao Yin, Paul Sajda and Stanley Osher
%%%              Columbia   , Rice,    , Columbia,      UCLA
%=========================================================================%

%%% 06/28/2007   Jianing Shi
%%% Logistic loss function

%%% 07/20/2007   Jianing Shi
%%% Modify the logistic loss function slightly
%%% Avoid numerical error due to floating point saturation

function y = LogisticLoss(x)

y = zeros(size(x));
indexA = find(x > -710);
y(indexA) = log(1+exp(-x(indexA)));
indexB = find(x < -709);
y(indexB) = log(1+exp(x(indexB)))-x(indexB);

