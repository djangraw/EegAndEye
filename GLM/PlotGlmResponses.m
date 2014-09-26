function h = PlotGlmResponses(varargin)

% Plots a topomovie for each regressor response function.
%
% h = PlotGlmResponses(results,tBinSize,color_limits)
% h = PlotGlmResponses(regressor_events,responseFns,tResponse,chanlocs,tBinSize,color_limits)
%
% INPUTS:
% - results is a struct with fields regressor_events, responseFns,
% tResponse, and NewEEG (which contains channel locations).
% - tBinSize and color_limits are plotting options for the MakeTopoMovie
% function (see that code for details).
% - regressor_events, responseFns, tResponse are the inputs and outputs of
% RunEegGlm (see that code for details).
% - chanlocs is the vector of channel structs usually in the 'chanlocs'
% field of an EEGLAB data struct.
%
% OUTPUTS:
% - h is a 1xn vector of structs with axis handles for the various parts of
% each figure (see MakeTopoMovie for details).
%
% Created 1/4/12 by DJ.
% Updated 4/30/13 by DJ - responseFns and regressor_events in cells
% Updated 3/19/14 by DJ - tResponse in cells

% Parse inputs
if nargin<=3
    a = varargin{1};
    a = UpdateGlmResultsFormat(a);
    regressor_events = a.regressor_events{a.iLevel};
    responseFns = a.responseFns{a.iLevel};
    tResponse = a.tResponse{a.iLevel};
    chanlocs = a.EEG.chanlocs;
    if nargin<2
        tBinSize = []; 
    else
        tBinSize = varargin{2};
    end
    if nargin<3
        color_limits = [];
    else
        color_limits = varargin{3};
    end
else
    regressor_events = varargin{1};
    responseFns = varargin{2};
    tResponse = varargin{3};
    chanlocs = varargin{4};
    % Set defaults
    if nargin<5
        tBinSize = [];
    else
        tBinSize = varargin{5};
    end
    if nargin<6
        color_limits = [];
    else
        color_limits = varargin{6};
    end
end


% Plot
nRegressors = size(responseFns,3);
for i=1:nRegressors
    h(i) = MakeTopoMovie(responseFns(:,:,i),tResponse,chanlocs,tBinSize,color_limits); % plot each regressor as a topomovie
    MakeFigureTitle(regressor_events{i},0); % figure title bar w/o text on figure
end

