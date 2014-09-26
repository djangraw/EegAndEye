function RF = ConvertBetaToRF(beta,T)

% ConvertBetaToRF(beta,T);
%
% INPUTS:
% - beta is a Dxp matrix
% - T is a scalar
%
% OUTPUTS: 
% - RF is a DxTxM matrix, where M=p/T is the number of event types.
%
% Created 9/19/14 by DJ.

D = size(beta,1);
p = size(beta,2);
M = p/T;
if mod(M,1)~=0
    error('p = %d is not evenly divisible by T = %d.',p,T);
end

RF = zeros(D,T,M);
for j=1:D
    for i=1:p
        % Get inidces
        iRF = ceil(i/T); % which contrast fn?
        iT = i-T*(iRF-1); % which time point in that CF?

        RF(j,iT,iRF) = beta(j,i);
    end    
end