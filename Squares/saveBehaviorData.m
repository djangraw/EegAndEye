function y = saveBehaviorData(subject,sessions,prefix)

% y = saveBehaviorData(subject,sessions)
%
% INPUTS:
% -subject is a scalar indicating the subject number.
% -sessions is an n-element vector indicating the session numbers to be loaded.
% -prefix is a string indicating the part of the filename before the first
% hyphen [default: 'sq']
%
% OUTPUTS:
% -y is an n-element vector of behavioral data structs.
% -a file called <prefix>-<subject>-all.mat will be saved in the current
% directory, and will contain the variable y.
%
% Created 12/13/11 by DJ.

if nargin<3
    prefix = 'sq';
end

nSessions = numel(sessions);

for i=1:nSessions
    load(sprintf('%s-%d-%d.mat',prefix,subject,sessions(i)));
    y(i) = x;
end

save(sprintf('%s-%d-all.mat',prefix,subject),'y');