function [RF_seq, t_seq] = ReconstructTrial(RF,tResponse,tStims,weights)

% [RF_seq, t_seq] = ReconstructTrial(RF,tResponse,tStims,weights)
%
% INPUTS:
% -RF is a DxTxP matrix in which row i, page j is the group response 
% function over time for event type j on electrode i.
% -tResponse is a T-element vector indicating the times (in ms) of each
% sample in the response functions.
% -tStims is an N-element vector indicating when each stimulus in the
% reconstructed sequence takes place.
% -weights is a n NxP matrix in which row k is the weights across RF
% that apply to stimulus k in the reconstructed sequence.
%
% OUTPUTS:
% -RF_seq is an DxM matrix of the response function in the reconstructed
% sequence.
% -t_seq is an M-element vector of the corresponding times (in ms).
%
% Created 8/6/14 by DJ.

% Declare constants
[D,T,P] = size(RF);
N = length(tStims);
% Get times in sequence
dt = tResponse(2)-tResponse(1); % assume constant dt
t_seq = (min(tStims) + min(tResponse)):dt:(max(tStims)+max(tResponse));
% Initialize sequence RF
M = length(t_seq);
RF_seq = zeros(D,M);
% Input weighted RFs into sequence
for i=1:N
    iStart = find(t_seq<=tStims(i),1,'last');
    iTimes = iStart-1+(1:T);
    
    for j=1:P
        RF_seq(:,iTimes) = RF_seq(:,iTimes) + weights(i,j)*RF(:,:,j);
    end
end