function [sigmasq_all, sigmasq_elec] = GetGlmMse(X,Y,beta)

% Calculate the mean squared error
%
% [sigmasq_all, sigmasq_elec] = GetMse(X,Y,beta)
%
% INPUTS:
% -X is nxp
% -Y is nxD
% -beta is Dxp
%
% OUTPUTS:
% -sigmasq_all is a scalar giving the MSE across all electrodes.
% -sigasq_elec is a 1xD vector of the MSE for each electrode.
%
% Created 8/27/14 by DJ.

% Get constants
[D,p] = size(beta);
n = size(X,1);

% Reconstruct data
reconstructed = (X*beta');
% Get residuals
resid = Y - reconstructed;

% Get MSE across all electrodes
sigmasq_all = sum(resid(:).^2)/((n-p)*D);

% Get MSE for each individual electrode
sigmasq_elec = nan(1,D);
for i=1:D
    sigmasq_elec(i) = sum(resid(:,i).^2)/(n-p);
end