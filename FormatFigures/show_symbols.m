function strOut = show_symbols(strIn)

% Edits the string (or cell array of strings) strIn so all its special 
% characters (especially underscores) will be displayed as-is on a plot, 
% not as subscripts, etc.
%
% strOut = show_us(strIn)
%
% INPUTS:
% - strIn is a string or cell array of string containing some characters
% that would usually be considered 'special' by the _printf functions.
%
% OUTPUTS:
% - strOut is the same as strIn, except escape characters are added so that
% no characters will be used as 'special'.
%
% Created 8/31/09 by DJ.
% Updated 9/1/09 by DJ.
% Updated 5/15/12 by DJ - added cell input support
% Updated 7/13/12 by DJ - comments

strOut = strIn;

if iscell(strIn) % for cell arrays of strings
    
    for i=1:numel(strIn) % recursively call function for each string
        strOut{i} = show_symbols(strIn{i});
    end
    
else % for strings (typical usage)
    
    invisible_chars = {'\','_','^'}; % SLASH MUST BE FIRST, IF IT'S USED!

    % Add slash before each invisible character
    for i=1:numel(invisible_chars)
        ind = strfind(strOut,invisible_chars{i});
        for j=numel(ind):-1:1
            strOut = [strOut(1:ind(j)-1) '\' strOut(ind(j):end)]; 
        end
    end
    
end