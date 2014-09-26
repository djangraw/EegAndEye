function [group_RF, group_p, group_SE] = RunTopLevelGlm_EEG(contrast_results,variances,multcorrect)

% [group_RF, group_p, group_SE] = RunTopLevelGlm_EEG(contrast_results,variances,multcorrect)
%  
% INPUTS:
% -contrast_results is a matrix of size (D,T,M,N), where D=#electrodes,
% T=#offsets, M=#contrasts, and N=#subjects. It contains the contrast
% values estimated by the level 1 GLM.
% -variances is a matrix of size (D,T,M,N). It contains the variance of the
% contrast values estimated by the level 1 GLM.
% - multcorrect is a string indicating the type of multiple-comparisons
% correction you want to use.  Supported: 'none' (default), 'fdr' (false
% discovery rate), 'bonferroni'.
%
% OUTPUTS:
% -group_RF is a matrix of size (D,T,M) indicating the values of the
% contrasts as learned by the GLM.
% -group_p is a matrix of the p values learned from the GLM (based on a t
% statistic, and NOT corrected for multiple comparisons).
% -group_SE is a matrix of the same size with the standard error of the
% group estimate.
%
% Created 4/18/13 by DJ.
% Updated 4/19/13 by DJ - added multcorrect input
% Updated 9/18/14 by DJ - added group_SE output

if nargin<3
    multcorrect = 'none';
end

% set up
[D,T,M,N] = size(contrast_results); % elecs, offsets, contrasts, subjects
p = T*M; % # "voxels"

% Set up
group_RF = nan(D,T,M);
group_p = nan(D,T,M);
group_SE = nan(D,T,M);

% Run GLMs
fprintf('---Beginning GLM on %d channels.\n',D);
parfor i=1:D
    tic;
    fprintf('channel %d/%d...',i,D);
    % reshape inputs
    betas = reshape(contrast_results(i,:,:,:),[p N]);
    sigma_sq = reshape(variances(i,:,:,:), [p N]);

    % Perform GLM
    [beta_g, p_g, se_g] = RunTopLevelGlm(betas, sigma_sq, multcorrect);

    % Reshape outputs
    group_RF(i,:,:) = reshape(beta_g,[1 T M]);
    group_p(i,:,:) = reshape(p_g,[1 T M]);
    group_SE(i,:,:) = reshape(se_g,[1 T M]);
    t = toc;
    fprintf('Done! took %.2f seconds.\n',t);
end
fprintf('--- Done!\n');

