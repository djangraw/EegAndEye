function data = SubtractOutComponents(data0, components)

% Find the projections of data onto given components and remove them.
%
% data = SubtractOutComponents(data0, components)
% 
% INPUTS:
% -data0 is a DxT matrix of data (channels x time).
% -components is a DxM vector of components to remove (chans x components). 
% Each column is a component (for example, HEOG or blink activity).
%
% OUTPUTS:
% -data is a DxT matrix (data0 with the given components projected out).
%
% Created 6/4/13 by DJ.

% declare shorter name for components input
A = components;

% Find pseudoinverse of components
V = inv(A'*A)*A'; 

% Find data with projection onto components removed
data = (eye(size(data0,1))-A*V)*data0;