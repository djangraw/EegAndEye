%%% 12/08/2010   Jianing Shi
%%% Parameter initialization for sparse logistic regression

function [Opts, Para] = SparseSolverInitialization()

%%% Choose classification evaluation
Opts.classopt = 3;

%%% Set parameters relating to sparse logistic regression
Para.maxStep = 100;
Para.lambda = 10;
Para.dt0 = 1e-3;