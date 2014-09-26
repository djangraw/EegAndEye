function y = loadBehaviorData(subject,sessions,prefix)

% y = loadBehaviorData(subject,sessions)
%
% INPUTS:
% -subject is a scalar indicating the subject number.
% -sessions is an n-element vector indicating the session numbers to be loaded.
% -prefix is a string indicating the part of the filename before the first
% hyphen [default: 'sq']
%
% OUTPUTS:
% -y is an n-element vector of behavioral data structs.
%
% Created 11/28/11 by DJ.
% Updated 12/4/11 by DJ - added prefix input
% Updated 3/6/13 by DJ - added || isempty(sessions)

if nargin<3
    prefix = 'sq';
end

if nargin==1 || isempty(sessions)
    load(sprintf('%s-%d-all.mat',prefix,subject)); % if it's saved
else
    nSessions = numel(sessions);

    for i=1:nSessions
        load(sprintf('%s-%d-%d.mat',prefix,subject,sessions(i)));
        y(i) = x;
    end
end