function [weights, courses, icawinv, varca] = ApplyIcaToGlmResults(h,tH,tRange,chanlocs,reg_names,nComps,colors)

% [weights, eigenvalues, icawinv, varca] = ApplyIcaToGlmResults(h,tH,tRange,chanlocs,reg_names,nComps,colors)
%
% INPUTS:
% -h is a DxTxN regressor, where D is the number of channels, T is the
%  number of time points, and N is the number of regressors.
% -tH is a T-element vector showing the corresponding time points relative
%  to a regressor event.
% -tRange is a 2-element vector indicating the start and end time of the
%  range you wish to use for calculating SVD. Default is [tH(1) tH(end)].
% -chanlocs (OPTIONAL) should be the chanlocs field of your EEGLAB data
%  struct, and should be included only if you wish to plot the top nComps
%  components found and their timecourses.
% -nComps is the number of components you want to include on the plot.
% -colors is an M-element cell array of strings or an Mx3 matrix indicating
%  the color for each line.
%
% OUTPUTS:
% -weights is a nComps x D matrix in which each row weights(i,:) is a set   
%  of electrode weights for component i.
% -courses is an NxTxD matrix in which courses(i,j,k) is the activity of
%  component k in response type i at time point j.
% -icawinv is a nCompsxD matrix in which each column icawinv(:,i) is a set of
%  inverse weights (for scalp map plotting) for component i.
% -varca is an nComps x 1 vector in which element i indicates the variance 
%  of component i's activity as fraction of the total data variance.
%
% Created 8/9/12 by DJ based on ApplySvdToGlmResults.
% Updated 2/15/13 by DJ - switched from RMS to RSS scaling, added subspace 
%   decomp to ica calls (project down to top nComps PC's, then find top 
%   nComps IC's in that subspace)
% Updated 3/13/14 by DJ - added colors input, icawinv, varca outputs.

disp('Finding components using ICA...')
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
if nargin<7
    colors = [];
end

% Get data in specified time range
D = size(h,1); % number of electrodes
N = size(h,3); % number of regressors
isIn = tH>=tRange(1) & tH<=tRange(2);
x = reshape(h(:,isIn,:),[D, sum(isIn)*N]); % data to use for ICA

avgref = 0;
if avgref
    disp('Re-referencing to average...');
    x = x-repmat(mean(x,1),size(x,1),1);
end
zeromean = 0;
if zeromean
    disp('Subtracting mean from ICA input...');
%     meanx = mean(x,2);
    x = x-repmat(mean(x,2),1,size(x,2));
%     h = h-repmat(meanx,[1,size(h,2),size(h,3)]);
end
% Convert to double for added precision
x = double(x);

% Perform subspace analysis
icatype = 'runica';
fprintf('Running %s\n...',icatype)
switch icatype
    case 'fastica'
        [icasig,icawinv,icaweights] = fastica(x,'lastEig',nComps,'numOfIC',nComps);
        icasphere = eye(size(icaweights,2));
    case 'runica'
        [icaweights,icasphere] = runica( x, 'lrate', 0.001, 'interupt', 'on', 'extended',1,'pca',nComps );
        icawinv = pinv(icaweights*icasphere);
end

scalerss = 1;
if scalerss
    disp('Scaling components to RSS microvolt');
    scaling = repmat(sqrt(sum(icawinv(:,:).^2))', [1 size(icaweights,2)]);
    icaweights = icaweights .* scaling;
    icawinv = pinv(icaweights * icasphere);
end

% Get variance of component activity in ICA window
varx = sum(var(x,0,2));
ca = icaweights*icasphere*x; % component activity
varca = var(ca,0,2)/varx; % variance of component activity as fraction of total data variance

sortcomps = 1;
if sortcomps
    disp('Resorting by variance of component activity in ica training window')    
    [varca,order] = sort(varca,'descend'); % put highest variance first
    icaweights = icaweights(order,:); % reorder weights
    icawinv = icawinv(:,order); % reorder winv's
end

weights = (icaweights*icasphere);
% Get activations
% icaact = (icaweights*icasphere)*x;

% Get timecourse of each component
courses = zeros(size(h,3),size(h,2),size(weights,1));
for i=1:size(weights,1)
    for j=1:size(h,3)
        courses(j,:,i) = weights(i,:)*h(:,:,j);
    end
end

if plotResults
    PlotSvdWeightsAndCourses(h,tH,weights,varca,chanlocs,reg_names,nComps,icawinv,colors);
end

disp('Done!')