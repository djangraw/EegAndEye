function h = DrawArc(P0,P1,P2,n)

% Draw an arc from P1 to P2 with its center at P0.
% 
% h = DrawArc(P0,P1,P2,n)
%
% INPUTS:
% -P0 is a 2x1 vector giving the x;y coordinate of the center of the arc.
% -P1 and P2 are 2x1 vectors giving the x;y coordinates of the arc's start
% and end.
% -n is the number of points in the arc [default=1000].
%
% OUTPUTS:
% -h is the handle of the plotted arc.
%
% (code adapted from
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/278048)
%
% Created 11/24/15 by DJ based on MATLAB thread 278048.


% Set up
if ~exist('n','var') || isempty(n)
    n = 1000; % The number of points in the arc
end
% Calculate
x0 = P0(1);
y0 = P0(2);
v1 = P1-P0;
v2 = P2-P0;
c = det([v1,v2]); % "cross product" of v1 and v2
a = linspace(0,atan2(abs(c),dot(v1,v2)),n); % Angle range
if c==0 % semicircle
    v3=[0 -1;1 0]*v1; % v1 rotated 90 degrees CCW
else
    v3 = [0,-c;c,0]*v1; % v3 lies in plane of v1 and v2 and is orthog. to v1
end
v = v1*cos(a)+((norm(v1)/norm(v3))*v3)*sin(a); % Arc, center at (0,0)
% Draw
h = plot(v(1,:)+x0,v(2,:)+y0); % Plot arc, centered at P0