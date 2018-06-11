function [struct1_crop,struct2_crop] = CropStructsToMatch(struct1,struct2,matchField)

% [struct1_crop,struct2_crop] = CropStructsToMatch(struct1,struct2,matchField)
%
% INPUTS:
% -struct1 and struct2 are structs with one field (matchField) in common.
% -matchField is a string indicating a field shared to struct1 and struct2.
% 
% OUTPUTS:
% -struct1_crop and struct2_crop are the input structs, but one has been
% cropped to match the other.
% 
% Created 7/26/16 by DJ.

if numel(struct1.(matchField))>numel(struct2.(matchField))
    if isequal(struct1.(matchField)(1:numel(struct2.(matchField))),struct2.(matchField))
        iMatchIn1 = 1:numel(struct2.(matchField));
    elseif isequal(struct1.(matchField)((end-numel(struct2.(matchField))+1):end),struct2.(matchField))
        iMatchIn1 = (1:numel(struct2.(matchField))) + numel(struct1.(matchField)) - numel(struct2.(matchField));
    else
        error('Structs could not be cropped to match!')
    end
    disp('Cropping struct1 to match struct2...')
    struct1_crop = CropStruct(struct1,iMatchIn1);
    struct2_crop = struct2;
elseif numel(struct2.(matchField))>numel(struct1.(matchField))
    if isequal(struct2.(matchField)(1:numel(struct1.(matchField))),struct1.(matchField))
        iMatchIn2 = 1:numel(struct1.(matchField));
    elseif isequal(struct2.(matchField)((end-numel(struct1.(matchField))+1):end),struct1.(matchField))
        iMatchIn2 = (1:numel(struct1.(matchField))) + numel(struct2.(matchField)) - numel(struct1.(matchField));
    else
        error('Structs could not be cropped to match!')
    end
    disp('Cropping struct2 to match struct1...')
    struct2_crop = CropStruct(struct2,iMatchIn2);
    struct1_crop = struct1;
else
    struct1_crop = struct1;
    struct2_crop = struct2;
end

% Check for success
if isequal(struct1_crop.(matchField), struct2_crop.(matchField))
    disp('Success!')
else
    error('Structs could not be cropped to match!')
end