function [erpsub time] = MakeErpMovie(ALLEEG,datadd,datsub,tbin_size,color_limits)

% Uses MakeTopoMovie specifically for a difference ERP from pop_comperp.
%
% [erpsub time] = MakeErpMovie( ALLEEG, datadd, datsub, tbin_size, color_limits)
% 
% Inputs (this list taken from pop_comperp's help file):
%   ALLEEG  - Array of loaded EEGLAB EEG structure datasets
%   datadd  - [integer array] List of ALLEEG dataset indices to average to make 
%             an ERP grand average and optionally to compare with 'datsub' datasets.
%   datsub  - [integer array] List of ALLEEG dataset indices to average and then 
%             subtract from the 'datadd' result to make an ERP grand mean difference. 
%             Together, 'datadd' and 'datsub' may be used to plot and compare grand mean 
%             responses across subjects or conditions. Both arrays must contain the same 
%             number of dataset indices and entries must be matched pairwise (Ex:
%             'datadd' indexes condition A datasets from subjects 1:n, and 'datsub', 
%             condition B datasets from the same subjects 1:n). {default: []}
%
% Optional Inputs (this list taken from MakeTopoMovie's help file):
% - tbin_size is the size of time bins (in seconds). Default is no
% binning in time.
% - color_limits is a 2-element vector indicating the min and max values
% you want to define the colorbar scale.  Default is 
% [max(data(:)) -max(data(:))].
% 
% Outputs (this list taken from pop_comperp's help file):
%   erpsub - Grand average (or rms) 'datadd' minus 'datsub' difference
%   times  - Vector of epoch time indices
%
% See pop_comperp for more info on how the ERP is calculated.
% See MakeTopoMovie for more info on how the data is plotted.
%
% Created 8/10/10 by DJ.
% Updated 11/2/10 by DJ - made positive voltages upward
% Updated 2/24/12 by DJ - added tbin_size, color_limits optional inputs
% Updated 3/1/12 by DJ - allow datasub to be empty

% handle defaults
if nargin<4
    tbin_size = [];
end
if nargin<5
    color_limits = [];
end

% Get and plot ERP
[erp1 erp2 erpsub time] = pop_comperp( ALLEEG, 1,datadd,datsub,'addavg','on','subavg','on','tplotopt',{'ydir',1});

% Decide which data to use
if isempty(erpsub)
    erpdata = erp1;
else
    erpdata = erpsub;
end

% Make TopoMovie figure
MakeTopoMovie(erpdata,time,ALLEEG(datadd(1)).chanlocs,tbin_size,color_limits);

