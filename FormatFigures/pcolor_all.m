function h = pcolor_all(varargin)

% Create a pseudocolor plot with pcolor that includes the outer edges of c
% and centers the boxes on the given x and y values.
%
% pcolor_all(x,y,c)
% pcolor_all(c)
%
% INPUTS:
% -c is a matrix of size mxn.
% -x is a monotonically increasing vector of length n.
% -y is a monotonically increasing vector of length m.
%
% OUTPUTS:
% -h is the handle of the resulting plot's surface object.
%
% Created 8/15/11 by DJ

% handle inputs
if nargin==1
    c = varargin{1};    
    x = 1:size(c,2);
    y = 1:size(c,1);
elseif nargin==2
    error(id('InvalidNumberOfInputs'),...
        'Must have one or three input data arguments.')
else
    x = varargin{1};
    y = varargin{2};
    c = varargin{3};
end

% make larger canvas
cnew = nan(size(c)+1);
cnew(1:end-1,1:end-1) = c;
% expand axes
xnew = interp1(1:numel(x),x,0.5:numel(x)+0.5,'linear','extrap');
ynew = interp1(1:numel(y),y,0.5:numel(y)+0.5,'linear','extrap');

% send results to pcolor
hh = pcolor(xnew,ynew,cnew);

if nargout==1
    h = hh;
end
