function [group_RF, group_p] = RunTopLevelGlm_EEG_flex(contrast_results,variances,group_x,group_contrast,multcorrect)

% RunTopLevelGlm_EEG_flex(contrast_results,variances,group_contrasts,multcorrect)
%  
% INPUTS:
% - contrast_results is a matrix of size (D,T,N), where D=#electrodes,
% T=#offsets, and N=#contrasts*subjects. It contains the contrast
% values estimated by the level 1 GLM.
% - variances is a matrix of size (D,T,N). It contains the variance of the
% contrast values estimated by the level 1 GLM.
% - group_x is an N x M matrix indicating the contribution of each
% subject/contrast to each of M group-level contributors. 
% - group_contrast is an M x Q matrix indicating which contrasts you
% want to compare.
% - multcorrect is a string indicating the type of multiple-comparisons
% correction you want to use.  Supported: 'none' (default), 'fdr' (false
% discovery rate), 'bonferroni'.
%
% OUTPUTS:
% - group_RF is a matrix of size (D,T,Q) indicating the values of the
% contrasts as learned by the GLM.
% - group_p is a matrix of the p values learned from the GLM (based on a t
% statistic, and corrected for multiple comparisons as indicated by the
% multcorrect input).
%
% Created 3/24/14 by DJ based on RunTopLevelGlm_EEG.

% Declare constants and default inputs
[D,T,N] = size(contrast_results); % elecs, offsets, level1contrasts*subjects
if nargin<3 || isempty(group_x)
    group_x = ones(N,1);
end
M = size(group_x,2);
if nargin<4 || isempty(group_contrast)
    group_contrast = eye(M); % one contrast for each group-level event type
end
Q = size(group_contrast,2);
if nargin<5 || isempty(multcorrect)
    multcorrect = 'none';
end

% Set up
% group_RF = nan(D,T,Q);
% group_p = nan(D,T,Q);
% Reshape & slice inputs
[betas,sigma_sq,group_RF_cell,group_p_cell] = deal(cell(1,D));
for i=1:D
    betas{i} = squeeze(contrast_results(i,:,:));
    sigma_sq{i} = squeeze(variances(i,:,:));
end


% Run GLMs
fprintf('---Beginning GLM on %d channels.\n',D);
parfor i=1:D
    tic;
    fprintf('channel %d/%d...',i,D);
    % reshape inputs
%     betas = squeeze(contrast_results(i,:,:)); %TxN matrix
%     sigma_sq = squeeze(variances(i,:,:)); %TxN matrix

    % Perform GLM
    [beta_g, p_g] = RunTopLevelGlm_flex(betas{i}, sigma_sq{i}, group_x, group_contrast, multcorrect);

    % Reshape outputs
    group_RF_cell{i} = reshape(beta_g,[1 T Q]);
    group_p_cell{i} = reshape(p_g,[1 T Q]);
    t = toc;
    fprintf('Done! took %.2f seconds.\n',t);
end
% Append all outputs
group_RF = cat(1,group_RF_cell{:});
group_p = cat(1,group_p_cell{:});

fprintf('--- Done!\n');

