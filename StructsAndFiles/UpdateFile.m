function UpdateFile(filename,variable,value)

% Updates the given variables in a file, then re-saves the file.
%
% UpdateFile(filename,variable,value)
% UpdateFile(filename,replacements)
%
% INPUTS:
% - filename is a string indicating the file you want to update.
% - variable is a string or cell array of strings indicating the variables
% to be updated.
% - value is the value you'd like to replace 'variable' with in the file.
% If 'variable' is a cell array, 'value' should be a cell array of the same
% size.
% -replacements (optional) is a struct in which each fieldname is a
% variable whose value will be set to the value of that field.
%
% Created 8/28/14 by DJ.

% Load into struct
fprintf('Loading %s...\n',filename);
R = load(filename);

% Update variables
fprintf('Updating variables...\n');
if ischar(variable)
    R.(variable) = value;
elseif iscell(variable)
    for i=1:length(variable)
        R.(variable{i}) = value{i};
    end
elseif isstruct(variable)
    fields = fieldnames(variable);
    for i=1:length(fields)
        R.(fields{i}) = variable.(fields{i});
    end
end

% Unpack struct into this workspace
Rfields = fieldnames(R);
UnpackStruct(R);

fprintf('Re-saving %s...\n',filename);
% Save results
save(filename,Rfields{:});