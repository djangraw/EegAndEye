function iBlinkSaccade = FindBlinkSaccades(blink_times, saccade_start_times, saccade_end_times)

% Gets the indices of saccades that are around blink events.
%
% iBlinkSaccade = FindHeogFixations(blink_times, saccade_start_times, saccade_end_times)
%
% INPUTS:
% -blink_times is a vector of times when blinks occurred.
% -saccade_start_times is an n-element vector of times when saccades began.
% -saccade_end_times is an n-element vector of times when saccades ended.
% saccade_start_times(i) and saccade_end_times(i) should refer to the start
% and end time of the same saccade.
%
% OUTPUTS:
% -iBlinkSaccade is a vector of indices of saccades in saccade_end_times
% that are around blink events.
%
% Created 6/4/13 by DJ.

% Set up
nBlinks = length(blink_times);
iBlinkSaccade = nan(1,nBlinks);

% Each blink event is inside a saccade, so...
% Find first saccade end time after each blink
for i=1:nBlinks
    iNextSac = find(saccade_start_times <= blink_times(i) & saccade_end_times >= blink_times(i),1);
    if ~isempty(iNextSac)
        iBlinkSaccade(i) = iNextSac;
    end
end

% Remove nan entries
iBlinkSaccade = iBlinkSaccade(~isnan(iBlinkSaccade));