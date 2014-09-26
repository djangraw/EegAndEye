function cmap = ThresholdedCmap(clim,cthresh,startmap)

% Create a colormap with low values set to zero
%
% cmap = ThresholdedCmap(clim,cthresh,startmap)
%
% INPUTS:
% - clim is a 2-element vector
% - cthresh is a scalar indicating the threshold below which colors should
% be set to the zero color.
% - startmap is the colormap you want to start with (default is jet(128))
%
% OUTPUTS:
% - cmap is the thresholded colormap, which can be used as input to
% colormap().
%
% Created 9/13/13 by DJ.

% Handle inputs
if nargin<3 || isempty(startmap)
    n = 128;
    startmap = jet(n);
elseif numel(startmap)==1
    n = startmap;
    startmap = jet(n);
else
    n = size(cmap,1);
end
cmap = startmap;

% Get values mapped to each color
cval = linspace(clim(1),clim(2),n); % value associated with each row of cmap

% Get zero and low values
i0 = find(cval>=0,1); % index of zero value
iLow = find(abs(cval)<cthresh); % indices below threshold

% Set new colormap
cmap(iLow, :) = repmat(cmap(i0,:),length(iLow),1); % set low values to 0 color