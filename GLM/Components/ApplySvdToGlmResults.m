function [weights, eigenvalues] = ApplySvdToGlmResults(h,tH,tRange,chanlocs,reg_names,nComps)

% [weights, eigenvalues] = ApplySvdToGlmResults(h,tH,tRange,chanlocs,reg_names,nComps)
%
% INPUTS:
% -h is a DxTxN regressor, where D is the number of channels, T is the
% number of time points, and N is the number of regressors.
% -tH is a T-element vector showing the corresponding time points relative
% to a regressor event.
% -tRange is a 2-element vector indicating the start and end time of the
% range you wish to use for calculating SVD. Default is [tH(1) tH(end)].
% -chanlocs (OPTIONAL) should be the chanlocs field of your EEGLAB data
% struct, and should be included only if you wish to plot the top nComps
% components found and their timecourses.
% -nComps is the number of components you want to include on the plot.
%
% OUTPUTS:
% -weights is a DxD matrix in which each column weights(:,i) is a set of  
% electrode weights for component i (NOTE: ica would call these 'icawinv')
% -eigenvalues is a D-element vector in which each element eigenvalues(i) 
% is the normalized eigenvalue for weights(:,i).
%
% Created 12/23/11 by DJ.
% Updated 1/6/12 by DJ - supports pooling over multiple regressors.
% Updated 7/18/12 by DJ - added nComps input

disp('Finding components using SVD...')
% Declare defaults
if nargin<2 || isempty(tH)
    tH = 1:size(h,2);
end
if nargin<3 || isempty(tRange)
    tRange = [tH(1), tH(end)];
end
if nargin<4 || isempty(chanlocs)
    plotResults = 0; 
else
    plotResults = 1; % plot the top components in the current figure
end
if nargin<5 
    reg_names = {};
end
if nargin<6
    nComps = 6;
end

% Get data in specified time range
D = size(h,1); % number of electrodes
N = size(h,3); % number of regressors
isIn = tH>=tRange(1) & tH<=tRange(2);
x = reshape(h(:,isIn,:),[D, sum(isIn)*N])'; % data to use for SVD

avgref = 0;
if avgref
    x = x-repmat(mean(x,2),1,size(x,2));
end

% Perform subspace analysis
[U,S,V] = svd(x,0);

% Declare outputs
weights = V;
eigenvalues = diag(S).^2;
eigenvalues = eigenvalues/sum(eigenvalues); % normalize so that each eigenvalue is the pct. variance explained by this component

if plotResults
    PlotSvdWeightsAndCourses(h,tH,weights',eigenvalues,chanlocs,reg_names,nComps,weights);
end

disp('Done!')