function [result, place] = find_first_word(text_file,word,format,stop_point)

% result = find_first_word(text_file,word,format,stop_line)
%
% This function reads a text file and finds all instances of a word.
% It returns a vector that has all of the strings after the word.
% Specifically designed for e-prime data text files.
% if option = 0, look through all file
% if option = 1 then look after practice trials (distinguished by 
% the first instance of "PracticeBlockList")
% format specifies the format of the output word_vec.
%   '%s' specifies a string: word_vec will be a cell array of strings.
%   '%d' specifies a double: word_vec will be a vector of numbers.
%
% Created 7/27/10 by DJ (based on find_word)
% Updated 8/31/10 by DJ - fixed numeric format bug.
% Updated 9/21/10 by DJ - added input stop_point and output place.

% Set defaults
if nargin<3
    format = '%s';
end

% Setup
fid = fopen(text_file);
fseek(fid,0,'eof'); % find end of file
eof = ftell(fid);
if nargin<4 || stop_point> (eof-length(word))
    stop_point = eof - length(word);  % using fscanf or fgetl after this point will cause an error
end
fseek(fid,0,'bof'); % rewind to beginning
%fseek(fid,400,-1);  % jump past irrelevant intro material (FOR EPRIME FILES!)

word_loc = '';  % current word/line being searched

% Main Loop
while ~strcmp(word_loc,word) && ftell(fid) < stop_point % if we haven't found it
  word_loc = fscanf(fid,'%s',1); % search next word
end
if strcmp(word_loc,word)  % if we've actually found the word 
    place = ftell(fid);
    result = fscanf(fid,format,1); % the next word will be the result we're looking for - read it in the proper format
else % if we didn't find the word
    result = [];
    place = [];
    warning(sprintf('Couldn''t find ''%s'' in document ''%s''!'...
        ,word,text_file))
end

% Clean up
fclose(fid); 

