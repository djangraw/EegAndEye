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
%%% Logistic loss function derivative

function y = LogisticLossDeri(x)

y  = -1./(exp(x)+1);