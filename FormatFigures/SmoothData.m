function data_out = SmoothData(data,sigma,halforfull)

% This function takes a vector of data as input and convolves it with a
% one-sided Gaussian kernel to smooth the curve.
%
% data_out = SmoothData(data,sigma,halforfull)
%
% INPUTS:
% - data: n-element vector of data to smooth
% - sigma: standard deviation of gaussian you wish to convolve with data
% - halforfull: string indicating whether a causal half-gaussian ('half')
% or a full gaussian ('full') should be used in the convolution.
%
% OUTPUTS:
% - data_out: n-element vector of data convolved with the half-gaussian
%
% Created 11/5/07 by DJ.
% Updated 11/6/07 by DJ.
% Updated 9/19/13 by DJ - interpolate, normalize, comments
% Updated 1/28/14 by DJ - fixed nan bug
% Updated 4/23/14 by DJ - added halforfull input

if nargin<2 || isempty(sigma)
    sigma = 1;
% 	sigma=numel(data)/10;    
end
if nargin<3 || isempty(halforfull)
    halforfull = 'half';
end

if sigma==0 % no smoothing
    data_out = data;
    return;
end

% call recursively for >2D data
if ndims(data)>2    
    data_out = nan(size(data));
    for i=1:numel(data(1,1,:));
        data_out(:,:,i) = SmoothData(data(:,:,i),sigma,halforfull);
    end
    return;
end

% The data we're smoothing may have NaNs involved.  
% Interpolate them for now, and we'll put the NaNs back later.
nan_indices = find(isnan(data));
num_indices = find(~isnan(data));
data(nan_indices) = interp1(num_indices,data(num_indices),nan_indices);

% Create narrow Gaussian kernel for smoothing the data
switch halforfull
    case 'half'
        dist = (0:(sigma*7))';
        gaus = [zeros(sigma*7,1); 1/(sigma*sqrt(2*pi))*exp(-((dist.^2)/(2*sigma^2)))]';
    case 'full'
        dist = (-round(sigma*7):round(sigma*7))';
        gaus = [1/(sigma*sqrt(2*pi))*exp(-((dist.^2)/(2*sigma^2)))]';
end    
gaus = gaus/sum(gaus); % normalize so that gaus sums to 1

% Smooth the data
newdata = conv2(data,gaus,'same');

% normalize the data
% newdata = newdata/max(newdata);

data_out = nan(size(data));
% data_out(nan_indices) = NaN;
data_out(num_indices) = newdata(num_indices);