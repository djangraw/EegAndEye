function b = GetContrastsFromResponseFns(a,option)

% Add and subtract various GLM response functions to get contrasts.
%
% b = GetContrastsFromResponseFns(a,option)
%
% INPUTS:
% - a is a struct containing fields 'regressor events' and 'responseFns',
% usually a GLM analysis you have saved and reloaded into a struct.
% - option is a string indicating what kind of contrasts you want
% (currently only '6con' is supported).
%
% OUTPUTS:
% - b is a copy of a with contrasts taking the place of regressors.
%
% Created 1/13/12

% Declare defaults
if nargin<2 || isempty(option)
    option = '6con';
end

% Set up
switch option
    case '6con'       
        b = a;
        b.regressor_events = {'Prep','Antic','Dist','Targ','Int','Com'};        
end
b.responseFns = zeros(size(a.responseFns,1),size(a.responseFns,2),numel(b.regressor_events));

% Separate out regressors
D0T = a.responseFns(:,:,strcmp('Dist-0T',a.regressor_events));
Integ = a.responseFns(:,:,strcmp('Integ',a.regressor_events));
D1T = a.responseFns(:,:,strcmp('Dist-1T',a.regressor_events));
Compl = a.responseFns(:,:,strcmp('Compl',a.regressor_events));

% Compute contrasts
for i=1:numel(b.regressor_events)
    switch b.regressor_events{i}
        case 'Prep'
            b.responseFns(:,:,i) = D0T + Integ - D1T - Compl;
        case 'Antic'
            b.responseFns(:,:,i) = -D0T - Integ + D1T + Compl;                
        case 'Dist'
            b.responseFns(:,:,i) = D0T - Integ + D1T - Compl;
        case 'Targ'
            b.responseFns(:,:,i) = -D0T + Integ - D1T + Compl;
        case 'Int'
            b.responseFns(:,:,i) = -D0T + Integ;
        case 'Com'
            b.responseFns(:,:,i) = -D1T + Compl;
    end
end