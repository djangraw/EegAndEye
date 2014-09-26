function newData = InterpolateBlinks(data,times,blinkRange)

% Linearly interpolate eye position or pupil size data that took place 
% during blinks
%
% InterpolateBlinks(data,times,blinkRange)
% InterpolateBlinks(data,times,x)
%
% INPUTS:
% -data is an n-element vector of eye position or pupil data.
% -times is an n-element vector of the corresponding times at which these
% data were taken (in eyelink samples, as in x.eyelink)
% -blinkRange is an mx2 matrix of blink start and end times.
% -x is a 3DSearch data struct.
%
% OUTPUTS:
% -newData is an n-element vector that is equal to data outside of blink
% ranges and interpolates the data inside the blink ranges.
%
% Created 8/8/11 by DJ.

% Handle inputs
if isempty(times)
    times = 1:length(data);
end
if isstruct(blinkRange)
    blinkRange = BlinkRange(blinkRange);
end
newData = data;
% If data is more than one column, have this function call itself recursively
if size(data,2)>1
    for i=1:size(data,2) % for each column in data
        newData(:,i) = InterpolateBlinks(data(:,i),times,blinkRange);
    end
    return;
end

% Set up
nBlinks = size(blinkRange,1);

% Interpolate using interp1
for i=1:nBlinks
    iStart = find(times<blinkRange(i,1),1,'last');
    iEnd = find(times>blinkRange(i,2),1,'first');
    newData(iStart:iEnd) = interp1([iStart iEnd], [data(iStart) data(iEnd)],iStart:iEnd);
end