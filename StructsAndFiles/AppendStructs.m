function combostruct = AppendStructs(structs_cell,cat_dim,fields_tocopy)

% Travels down struct fields until reaching non-struct value fields, then 
% appends values from individual structs to create one 'combo-struct'.
%
% combostruct = AppendStructs(structs_cell,cat_dim,fields_tocopy)
%
% INPUTS:
% -structs_cell is a cell array or vector of structs that need not share
% all fields in common.
% -cat_dim is the dimension in which numeric fields should be concatenated.
% For example, cat_dim=1 indicates concatenation in rows.
% - fields_tocopy is a cell array of strings indicating which fields should
% be included in the final combo-struct. All sub-fields of these fields 
% will be copied. [default: union of all fields in all structs].
%
% OUTPUTS:
% -combostruct is a struct with fields listed in fields_tocopy.
%
% Created 9/22/14 by DJ.
% Updated 9/23/14 by DJ - comments.

% Set up
if ~exist('cat_dim','var') || isempty(cat_dim)
    cat_dim = 1; % append rows
end
if ~exist('fields_tocopy','var')
    fields_tocopy = [];
end
nStructs = numel(structs_cell);

% Convert struct vector to cell array for maximum flexibility
if ~iscell(structs_cell)
    structs_vec = structs_cell;
    structs_cell = cell(1,nStructs);
    for iStruct=1:nStructs
        structs_cell{iStruct} = structs_vec(iStruct);
    end
    clear structs_vec
end

% Automatically determine fields to copy (if not specified)
if isempty(fields_tocopy)
    % Get a list of all the fields
    structfields = cell(1,nStructs);
    for iStruct=1:nStructs
        structfields{iStruct} = fieldnames(structs_cell{iStruct});    
    end
    % Find a unique list of all the fields
    fields_tocopy = unique([structfields{:}]);
end

% Fill in each field
nFields = numel(fields_tocopy);

for iField = 1:nFields    
    % initialize
    thisfield = fields_tocopy{iField};
    combostruct.(thisfield) = [];
    % Travel down structs to get to numbers or cells
    if isstruct(structs_cell{1}.(thisfield))
        substructs_cell = cell(1,nStructs);
        for iStruct = 1:nStructs
            if isfield(structs_cell{iStruct},thisfield)
                substructs_cell{iStruct} = structs_cell{iStruct}.(thisfield);
            else
                substructs_cell{iStruct} = struct();
            end
        end  
        % Call recursively
        combostruct.(thisfield) = AppendStructs(substructs_cell,cat_dim);
    else
        % Append numbers/strings/cells
        for iStruct = 1:nStructs     
            if isfield(structs_cell{iStruct},thisfield)
                combostruct.(thisfield) = cat(cat_dim, combostruct.(thisfield), structs_cell{iStruct}.(thisfield));
            end
        end
    end
end

