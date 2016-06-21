function outputStruct = ApplyDefaultsToStruct(inputStruct,defaultStruct)

% outputStruct = ApplyDefaultsToStruct(inputStruct,defaultStruct)
%
% INPUTS:
% -inputStruct is a struct with some, but not all, of the fields the
% calling function requires.
% -defaultStruct is a struct filled with all the fields the calling
% function requires.
% 
% OUTPUT:
% -outputStruct is a struct with all the fields found in either input. If
% both inputStruct and defaultStruct have a field, the value in inputStruct
% will be used. Otherwise, the value in defaultStruct will be used.
%
% Created 6/16/16 by DJ.

% Handle empty input
if isempty(inputStruct)
    outputStruct = defaultStruct;
    return;
end
% Set up
outputStruct = inputStruct;
fields = fieldnames(defaultStruct);

% Loop through each field
for i=1:numel(fields)
    if ~isfield(inputStruct,fields{i})
        % Fill field with default value
        outputStruct.(fields{i}) = defaultStruct.(fields{i});
    elseif isstruct(inputStruct.(fields{i}))
        % Call recursively
        outputStruct.(fields{i}) = ApplyDefaultsToStruct(inputStruct.(fields{i}), defaultStruct.(fields{i}));
    end
end