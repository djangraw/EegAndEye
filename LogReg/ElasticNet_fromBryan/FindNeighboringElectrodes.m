function [closest, distance] = FindNeighboringElectrodes(EEG)

% For each electrode, sorts the other electrodes by distance from it.
%
% [closest, distance] = FindNeighboringElectrodes(EEG);
%
% Also displays the 4 closest electrodes in a format consistent with
% Bryan's function computeElectrodeTptLaplacianMatrix, to be pasted in if
% you have a new cap that's unfamiliar to that function.
%
% INPUTS:
% -EEG is a standard EEGLAB data structure with the field 'chanlocs'
% initialized (i.e. a .loc file has been imported).
% 
% OUTPUTS:
% -closest is an nxn matrix of electrode indices (where n is the # of 
% electrodes). Each row i lists the indices of the electrodes, in the order
% of how close it is.
% -distance is an nxn matrix containing the euclidean distances used to 
% sort the electrodes.  distance(i,j) is the distance between electrode i 
% and electrode closest(i,j).
%
% Created 6/2/11 by DJ.

% Get positions of electrodes in Cartesian coordinates
x = [EEG.chanlocs(:).X];
y = [EEG.chanlocs(:).Y];
z = [EEG.chanlocs(:).Z];
% get # of electrodes
n = numel(EEG.chanlocs);
% set up
distance = zeros(n);
closest = zeros(n);

% Main loop
for i=1:n
% Get Euclidean distances
d = sqrt((x(i)-x).^2 + (y(i)-y).^2 + (z(i)-z).^2); % vector of distances from electrode i
[dsorted, order] = sort(d,'ascend'); % sort distances
% add to output variables
distance(i,:) = dsorted;
closest(i,:) = order;

% print in format of Bryan's functiopn computeElectrodeTptLaplacianMatrix
fprintf('AC(%d,[%d,%d,%d,%d]) = 1;\n',closest(i,1), closest(i,2), closest(i,3), closest(i,4), closest(i,5));

end