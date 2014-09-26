function [RF_seq,t_seq,titlestring] = PlotReconstructedTrial(RF,tResponse,event_list, event_seq,chanlocs,chansToPlot)

% [RF_seq,t_seq] = PlotReconstructedTrial(RF,tResponse,event_list, event_seq,chanlocs,chansToPlot)
%
% INPUTS:
% -RF is a DxTxP matrix in which row i, page j is the group response 
% function over time for event type j on electrode i.
% -tResponse is a T-element vector indicating the times (in ms) of each
% sample in the response functions.
% -event_list is a P-element cell array of strings indicating the names of
% the events in the RF matrix.
% -event_seq is an N-element cell array of strings describing the events in
% the sequence you'd like to reconstruct (usually 5 squares and a circle).
% -chanlocs is a D-element vector of structs
% -chansToPlot
%
% OUTPUTS:
% -RF_seq is an DxM matrix of the response function in the reconstructed
% sequence.
% -t_seq is an M-element vector of the corresponding times (in ms).
%
% Created 8/6/14 by DJ

[D,T,P] = size(RF);
N = length(event_seq);
tStims = 0:500:((N-1)*500);

rampweights = (1:5)/mean(1:5);

weights = zeros(N,P);
for i=1:N
    isPrimaryEvent = strcmp(event_seq{i},event_list);
    isRampEvent = strcmp([event_seq{i} '-RampUp'],event_list);
        
    weights(i,isPrimaryEvent) = 1;
    if i<=length(rampweights)
        weights(i,isRampEvent) = rampweights(i);
    end
end
% Add in TrialStart event
isTrialStartEvent = strcmp('TrialStart',event_list);
weights(1,isTrialStartEvent) = 1;

% Add in Square event
isSquareEvent_list = strcmp('sf-Square',event_list) | strcmp('Square',event_list) | strcmp('sf3-Square',event_list);
isSquareEvent_seq = ~(strcmp('sf-Circle',event_seq) | strcmp('Circle',event_seq) | strcmp('sf3-Circle',event_seq));
weights(isSquareEvent_seq,isSquareEvent_list) = 1;

[RF_seq, t_seq] = ReconstructTrial(RF,tResponse,tStims,weights);

% PLOT!
clf;
PlotResponseFnsGrid(RF_seq,{'sequence'},t_seq,chanlocs,chansToPlot);
titlestring = [];
for i=1:N 
    titlestring = [titlestring event_seq{i} ', '];
end
titlestring = titlestring(1:end-2);
% annotate
MakeFigureTitle(['Reconstructed Sequence:\n' titlestring]);