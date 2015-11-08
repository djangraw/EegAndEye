function croppedStruct = CropStruct(structToCrop,indicesToKeep,fieldsToCrop)

% CropStruct(structToCrop,indicesToKeep,fieldsToCrop)
%
% INPUTS:
% - structToCrop is a struct.
% - indicesToKeep is a vector of binary values of the same size as
% structToCrop.(fieldsToCrop{i}) OR a vector of indices whose max value is
% less than numel(structToCrop.(fieldsToCrop{i})).
% - fieldsToCrop is a cell array of strings indicating the fields of
% structToCrop that you'd like to crop. [default: all fields]
%
% OUTPUTS:
% -croppedStruct is a struct with the same fields as structToCrop.
%
% Created 9/14/15 by DJ.

% Handle defaults
if ~exist('fieldsToCrop','var') || isempty(fieldsToCrop)
    % get field names
    fieldsToCrop = fieldnames(structToCrop);
end

% First, copy all of structToCrop. Include the fields you don't want to crop.
croppedStruct = structToCrop;

% crop each field in turn
for i=1:numel(fieldsToCrop)    
    croppedStruct.(fieldsToCrop{i}) = structToCrop.(fieldsToCrop{i})(indicesToKeep);
end