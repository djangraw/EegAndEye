function [C_crop, tContrast] = BuildSequenceContrast(event_list,tResponse,event_seq,tEvents,tContrast)

% Created 8/11/14 by DJ.

dt = tResponse(2)-tResponse(1); % assume constant
tSequence = (min(tResponse)+min(tEvents)):dt:(max(tResponse)+max(tEvents));

if ~exist('tContrast','var')
    tContrast = tSequence;
end

T = length(tResponse);
P = length(event_list);
N = numel(tEvents);

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

% Initialize sequence RF
M = length(tSequence);
C_seq = zeros(P*T,M);
% Input weighted identity matrices into sequence
for i=1:N
    iStart = find(tSequence<=tEvents(i)+tResponse(1),1,'last');
    iTimes = iStart-1+(1:T);
    
    for j=1:P
        C_seq((j-1)*T+(1:T),iTimes) = C_seq((j-1)*T+(1:T),iTimes) + weights(i,j)*eye(T);
    end
end

% Crop
iStart = find(tSequence<=tContrast(1),1,'last');
iEnd = find(tSequence<=tContrast(end),1,'last');
iTimes = iStart:iEnd;
if numel(iTimes)~=numel(tContrast)
    error('Times in sequence do not match number of elements in tContrast!')
end
C_crop = C_seq(:,iTimes);




