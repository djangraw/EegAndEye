function subplots = GetSubplots(fignums)

% Returns the handles for all axes on the given figure(s).
%
% subplots = GetSubplots(fignums)
%
% INPUTS:
% - fignums is a vector of figure handles.
%
% OUTPUTS:
% - subplots is a column vector of the axis handles.
%
% Created 3/5/14 by DJ.

subplots = findall(fignums,'Type','axes','Tag','');