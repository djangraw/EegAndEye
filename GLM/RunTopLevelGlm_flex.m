function [beta_g, p_g] = RunTopLevelGlm_flex(level1beta, level1beta_var, group_x, group_contrast, multcorrect)

% [beta_g, p_g] = RunTopLevelGlm_flex(level1beta, level1beta_var, group_x, group_contrast, multcorrect)
%
% INPUTS:
% - level1beta is a T x N matrix
% - level1beta_var is a T x N matrix
% - group_x is an N x M matrix indicating the contribution of each
% subject/contrast to each of M group-level contributors. 
% - group_contrast is an M x Q matrix indicating which contrasts you
% want to compare.
% - multcorrect is a string indicating the type of multiple-comparisons
% correction you want to use.  Supported: 'none' (default), 'fdr' (false
% discovery rate), 'bonferroni'.

%
% OUTPUTS:
% - beta_g is a T x M vector of group-level betas.
% - p_g is a T x M vector of group-level p values based on the t statistic.
%
% Created 3/24/14 by DJ based on RunTopLevelGlm.

[T,N] = size(level1beta); % # "voxels", # subjects
if nargin<3 || isempty(group_x)
    group_x = ones(N,1);
end
M = size(group_x,2);
if nargin<4 || isempty(group_contrast)
    group_contrast = eye(M); % one contrast for each group-level event type
end
if nargin<5 || isempty(multcorrect)
    multcorrect = 'none';
end

% Set up
beta_star = reshape(level1beta,[T*N,1]); % output of model
sigma_sq = reshape(level1beta_var,[T*N,1]); % lower-level variances
% Construct group regressor matrix 
Xg = zeros(N*T,M*T); % group regressor matrix
for i=1:N
    for j=1:M
        Xg((i-1)*T+(1:T),(j-1)*T+(1:T)) = group_x(i,j)*eye(T);
    end
end
% Construct group contrast matrix
Q = size(group_contrast,2);
C_g = zeros(M*T,Q*T);
for j=1:M
    for k=1:Q
        C_g((j-1)*T+(1:T),(k-1)*T+(1:T)) = group_contrast(j,k)*eye(T);
    end
end
% C_g = eye(T); % simple contrast

% estimate group betas
beta_g_hat = (Xg'*Xg)^(-1) * Xg'*beta_star;

% get group level variance
% sigma_g_hat_sq = (Xg*beta_g_hat - beta_star)'*(Xg*beta_g_hat - beta_star);
dof = N*T - M*T;
sigma_g_hat_sq = (Xg*beta_g_hat - beta_star)'*(Xg*beta_g_hat - beta_star)/dof;
% sigma_g_hat_sq = sigma_g_hat_sq * c_g*(Xg'*Xg)^(-1)*c_g'; % NEW 11/11/13
Wg_vec = 1./sqrt(sigma_sq + sigma_g_hat_sq); % vector version of Wg
Wg = diag(Wg_vec); % weighting/whitening matrix

% re-estimate group betas
beta_g = (Xg'*Wg'*Wg*Xg)^(-1) * (Xg'*Wg'*Wg*beta_star);

% Get variance of contrasts
% resid_sq = (Wg*Xg*beta_g - Wg*beta_star)'*(Wg*Xg*beta_g - Wg*beta_star)/dof; % NEW 11/11/13
sigma_g_sq = nan(M*T,1);
for k=1:M*T
    sigma_g_sq(k) = C_g(k,:)*(Xg'*Wg'*Wg*Xg)^(-1)*C_g(k,:)';
end

% get std err of contrasts
% s_g = sqrt(sigma_g_sq/dof);
s_g = sqrt(sigma_g_sq);

% get t statistic
t_g = C_g*beta_g./s_g;

% get p values
Pval_start = tcdf(t_g,dof);
% Apply multiple-comparisons correction
switch multcorrect
    case 'none'
        Pval = Pval_start;
    case 'fdr' % false discovery rate
        isHigh = Pval_start>0.5;
        Pval_start(isHigh) = 1-Pval_start(isHigh);        
        Pval = reshape(mafdr(Pval_start(:),'bhfdr',true),size(Pval_start));
        Pval(isHigh) = 1-Pval(isHigh);
    case 'bonferroni' % bonferroni correction
        isHigh = Pval_start>0.5;
        Pval_start(isHigh) = 1-Pval_start(isHigh);
        Pval = bonf_holm(Pval_start);
        Pval(Pval>0.5) = 0.5;
        Pval(isHigh) = 1-Pval(isHigh);
    otherwise
        error('multiple comparisons method not recognized!');
end
p_g = Pval;