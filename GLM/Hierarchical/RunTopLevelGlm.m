function [beta_g, p_g, s_g] = RunTopLevelGlm(betas, beta_var, multcorrect)

% [beta_g, p_g, s_g] = RunTopLevelGlm(betas, beta_var, multcorrect)
%
% INPUTS:
% - betas is a p x N matrix
% - sigma_sq is a p x N matrix
% - multcorrect is a string indicating the type of multiple-comparisons
% correction you want to use.  Supported: 'none' (default), 'fdr' (false
% discovery rate), 'bonferroni'.
%
% OUTPUTS:
% - beta_g is a p x 1 vector of group-level betas.
% - p_g is a p x 1 vector of group-level p values based on the t statistic.
% - s_g is a p x 1 vector of group-level stderr values.
%
% Created 4/18/13 by DJ.
% Updated 4/19/13 by DJ - debugged, added multcorrect input
% Updated 8/5/14 by DJ - store XWinv so we don't have to recompute.
% Updated 9/18/14 by DJ - added std err output

if nargin<3 || isempty(multcorrect)
    multcorrect = 'none';
end

% Set up
[p,N] = size(betas); % # "voxels", # subjects
beta_star = reshape(betas,[p*N,1]); % output of model
Xg = repmat(eye(p),N,1); % group regressor matrix
sigma_sq = reshape(beta_var,[p*N,1]); % lower-level variances
c_g = eye(p); % simple contrast

% estimate group betas
beta_g_hat = (Xg'*Xg)^(-1) * Xg'*beta_star;

% get group level variance
% sigma_g_hat_sq = (Xg*beta_g_hat - beta_star)'*(Xg*beta_g_hat - beta_star);
dof = p*N-p;

% for i=1:10
%     fprintf('iter %d...',i);

    sigma_g_hat_sq = (Xg*beta_g_hat - beta_star)'*(Xg*beta_g_hat - beta_star)/dof;
    % sigma_g_hat_sq = sigma_g_hat_sq * c_g*(Xg'*Xg)^(-1)*c_g'; % NEW 11/11/13
    Wg_vec = 1./sqrt(sigma_sq + sigma_g_hat_sq); % vector version of Wg
    Wg = diag(Wg_vec); % weighting/whitening matrix

    % re-estimate group betas
    XWinv = (Xg'*Wg'*Wg*Xg)^(-1);
    beta_g = XWinv * (Xg'*Wg'*Wg*beta_star);

    % Get variance of contrasts
    % resid_sq = (Wg*Xg*beta_g - Wg*beta_star)'*(Wg*Xg*beta_g - Wg*beta_star)/dof; % NEW 11/11/13
    sigma_g_sq = nan(p,1);
    for k=1:p
        sigma_g_sq(k) = c_g(k,:)*XWinv*c_g(k,:)';
    end

    % get std err of contrasts
    % s_g = sqrt(sigma_g_sq/dof);
    s_g = sqrt(sigma_g_sq);
    
    % get t statistic
    t_g = c_g*beta_g./s_g;
    
%     beta_g_hat = beta_g;
%     % Print results
%     fprintf('Mean stderr = %g, t = %g\n',mean(s_g(:)),mean(t_g(:)));
% end

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