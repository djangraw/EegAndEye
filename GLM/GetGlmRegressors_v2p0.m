function [regressors, S] = GetGlmRegressors_v2p0(t,tEvents,tArtifacts,regressor_range,artifact_range,weights,stddev,zeroPadS)

% [regressors, S] = GetGlmRegressors_v2p0(t,tEvents,tArtifacts,regressor_range,weights,stddev,zeroPadS)
%
% INPUTS:
% - t is a vector of length L indicating the times at which EEG samples 
% were recorded.
% - tEvents is a cell array of length Nr, in which tEvents{i} is a vector 
% of all the times when event i happened.
% - tArtifacts is a vector of all the times when events that may cause an 
% artifact happened.  Times within artifact_range samples of these events 
% will set to zero in the regressors or S output matrix (see below).
% - regressor_range is a 2-element vector indicating how many EEG samples 
% before nd after each event will be affected by the event. For example, 
% [-2 5] indicates that 2 samples before to 5 samples after will be used. 
% [6indicates that 6 samples before to 6 samples after will be used.
% - artifact_range is a 2-element vector indicating how many EEG samples 
% before and after each artifact event should be zeroed out in the
% regressors or S matrix. Format is the same as regressor_range.
% - weights is a cell array of length Nr, in which weights{i} is a vector
% of the desired regressor values when event i happened. (default = ones)
% - stddev is a scalar indicating the standard deviation (in samples) of 
% the gaussian you wish to convolve with the results.  Default is 0, 
% indicating a delta function at the time of each event.
% - zeroPadS is a boolean value indicating whether the S matrix should be
% zero padded around the artifacts.  Otherwise, the regressors will be zero
% padded, eliminating any events that are too close to an artifact.
%
% OUTPUTS:
% - regressors is an Nr x L matrix in which each row i is 1 at times when
% event i happened.
% - S is the regressor matrix for input into SimpleGLM.m.  It will be of
% size L x ((diff(regressor_range)+1)*Nr).
%
% Created 12/22/11 by DJ (called GetGlemRegressors.m).
% Updated 6/18/12 by DJ - added weights input.
% Updated 6/26/12 by DJ - s,S will be booleans if weights are binary
% Updated 7/16/12 by DJ - added gaussian convolution option, stddev input
% Updated 3/22/13 by DJ - use asymmetric influence, separate artifact 
%  influence. (_v2p0)

% Handle inputs
if nargin<5 || isempty(artifact_range)
    warning('artifact_range not specified!');
    artifact_range = regressor_range;
end
if nargin<6 || isempty(weights) % default to all ones
    weights = cell(size(tEvents));
end
if nargin<7 || isempty(stddev)
    stddev = 0;
end
if nargin<8 || isempty(zeroPadS)
    zeroPadS = true;
end
% Check error conditions
if numel(regressor_range)>2 || numel(artifact_range)>2
    error('range inputs can only be length 1 or 2!');
end

% Handle 1-element ranges
if numel(regressor_range) == 1
    regressor_range = [-regressor_range, regressor_range];
end
if numel(artifact_range) == 1
    artifact_range = [-artifact_range, artifact_range];
end


% Declare constants
Nr = numel(tEvents);
L = length(t);
Nh = diff(regressor_range)+1; % length of response vectors

% Fill in any empty weights cells with default
for i=1:numel(weights)
    if isempty(weights{i})
        weights{i} = ones(size(tEvents{i}));
    end
end
% Check if weights are binary (so we can save memory)
if all([weights{:}]==1) && stddev==0
    useBinary = true;
else
    useBinary = false;
end


% Regressors
if useBinary
    regressors = false(Nr,L); 
else
    regressors = zeros(Nr,L); % each row is a regressor
end
iEvents = cell(size(tEvents));
if numel(unique(diff(t)))==1;
    disp('fast way')
    % find indices of regressor events
    dt = t(2)-t(1);
    for i=1:Nr                
        iEvents{i} = round(tEvents{i}/dt);
    end
    % Find indices of artifact events
    iArtifacts = round(tArtifacts/dt);
        
    
else
    disp('slow way')
    % find indices of regressor events
    for i=1:Nr
        iEvents{i} = nan(1,length(tEvents{i}));
        for j=1:length(tEvents{i})
            iEvents{i}(j) = find(t>=tEvents{i}(j),1);
        end
    end
    % Find indices of artifact events
    iArtifacts = nan(size(tArtifacts));
    for j=1:numel(tArtifacts)
        iArtifacts(j) = find(t>=tArtifacts(j),1);
    end
end

% Make regressor vector
for i=1:Nr
    regressors(i,iEvents{i}) = weights{i};
end

% Find times around artifact events
isBad = false(1,size(regressors,2));
for j=1:numel(iArtifacts)        
    if iArtifacts(j)<=-artifact_range(1)
        isBad(1:iArtifacts(j)+artifact_range(2)) = true;
    elseif iArtifacts(j)>L-artifact_range(2)
        isBad(iArtifacts(j)+artifact_range(1):end) = true;
    else
        isBad(iArtifacts(j)+(artifact_range(1):artifact_range(2))) = true;
    end
end

% Zero out regressors around artifact events
if ~zeroPadS
    regressors(:,isBad) = 0;
end


% Gaussian stamp
if stddev>0
    % Create gaussian stamp
    width = ceil(abs(norminv(0.01,0,stddev)));
    x = -width:width;
    y = normpdf(x,0,stddev); % do I need to use subtraction on normcdf values?
    z = y/sum(y); % normalize to sum to 1
%     plot(x,z)

    % Convolve regressors with gaussian stamp
    bigconv = conv2(regressors,z);
    regressors = bigconv(:,(1+(length(z)-1)/2) : (end-(length(z)-1)/2)); % get rid of padding
end

% free up memory
% clear iArt*
clear iEvents
clear t*
clear weights 
clear width x y z bigconv

% Construct regressor matrix
if useBinary
    S = false(L,Nh*Nr); 
else 
    S = zeros(L,Nh*Nr); % Changed 6/18/12 to make non-binary regressors work
end
for i=1:Nh
    ih = i+regressor_range(1)-1; % samples relative to actual event time
    if ih<0
        iS1 = 1;
        ir1 = 1-ih;
        iS2 = L+ih;
        ir2 = L;
    elseif ih>0
        iS1 = 1+ih;
        ir1 = 1;
        iS2 = L;
        ir2 = L-ih;
    else
        iS1 = 1;
        ir1 = 1;
        iS2 = L;
        ir2 = L;   
    end
    S(iS1:iS2,i:Nh:end) = regressors(:,ir1:ir2)';
end
S = sparse(S);

% Zero out S matrix around artifact events
if zeroPadS
    S(isBad,:) = 0;
end
    

disp('Done!')


