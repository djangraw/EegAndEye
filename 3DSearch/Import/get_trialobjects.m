function [trialobjects lifetimes] = get_trialobjects(text_file)

% trialobjects = get_trialobjects(text_file)
% This function reads an Eyelink text file from a 3D search (3DS) task and
% parses information about the objects created in that task into useful
% structs for later analysis.
% - Input text_file is the filename of the Unity Log (something like
% '3DS-0-0.txt').
% - Output trialobjects will be an n-element array of cells, where n is the
% number of trials in the file.  Each cell will be an m-element array of
% structs, with each struct trialobjects{i}(j) containing info about object
% #j of trial #i.
%
% Created 6/15/10 by DJ.
% Updated 6/24/10 by DJ - use Eyelink file instead of Unity Log.
% Updated 7/22/10 by DJ - formatting of textscan line has changed.
% Updated 7/26/10 by DJ - fixed strings in struct, handles Destroy All
% Updated 7/29/10 by DJ - store time info (system/Hz-dependent) separately
%   from object info (system-independent).
% Updated 8/25/10 by DJ - require only 'END TRIAL', not 'Destroyed All
%   Objects'.
% Updated 4/13/11 by DJ - added rotation information (started being logged
%   in 3DS version 5.8).
% Updated 10/26/11 by DJ - add support for non-Clone objects

% Setup
fid = fopen(text_file);
fseek(fid,0,'eof'); % find end of file
eof = ftell(fid);
fseek(fid,0,'bof'); % rewind to beginning

trialobjects = {}; %each cell is a vector of structs
lifetimes = {}; % each cell is an nObjects x 2 matrix
trial = 0;          % number of trials

% Find Trial markers and read in object info.
while ftell(fid) < eof % if we haven't reached the end of the text file
    % Find the start of a new trial
    str = fgetl(fid); % read in next line of text file
    if findstr(str,'LOAD TRIAL') % check for the code-word indicating loading started
        trial = trial+1; % We've reached the next trial, so increment the trial counter
        
        % Get all objects created for this trial
        while isempty(findstr(str,'END TRIAL')) % until the trial ends
            str = fgetl(fid); % read in next line of text file
            if findstr(str,'Created Object')
                if ~isempty(findstr(str,'Clone'))
                    objectinfo = textscan(str,'MSG %d Created Object # %d %s Clone %s %s %f %f %f',... % read info about object
                    'Delimiter',',() ','MultipleDelimsAsOne',true); % read info about object (line 2 of 2)
                        % Delimiter field makes us stop reading strings at those
                        % characters (important for the NAME(Clone) objects).
                else
                    objectinfo = textscan(str,'MSG %d Created Object # %d %s %s %s %f %f %f',... % read info about object
                    'Delimiter',',() ','MultipleDelimsAsOne',true); % read info about object (line 2 of 2)
                end
                    
                % Read the info we just parsed using sscanf into the proper variables
                objecttime = objectinfo{1}; 
                objectnum = objectinfo{2};
                objectname = objectinfo{3}{1};
                objecttag = objectinfo{4}{1};
                objecttype = objectinfo{5}{1};
                objectpos = [objectinfo{6:8}];
                if length(objectinfo)>8
                    objectrot = [objectinfo{9:12}];
                else
                    objectrot = [0 0 0 1];
                end
                lifetimes{trial}(objectnum,1:2) = [objecttime, 0]; % create time, destroy time
                trialobjects{trial}(objectnum).name = objectname;
                trialobjects{trial}(objectnum).type = objecttype;
                trialobjects{trial}(objectnum).tag = objecttag;
                trialobjects{trial}(objectnum).createposition = objectpos;
                trialobjects{trial}(objectnum).createrotation = objectrot;
            elseif findstr(str,'Destroyed Object')
                objectinfo = textscan(str,'MSG %d Destroyed Object # %d %s Clone %s %s',... % read info about object
                'Delimiter',',() ','MultipleDelimsAsOne',true); % read info about object (line 2 of 2)
                    % Delimiter field makes us stop reading strings at those
                    % characters (important for the NAME(Clone) objects).
                % Read the info we just parsed using sscanf into the proper variables
                objecttime = objectinfo{1}; 
                objectnum = objectinfo{2};
                lifetimes{trial}(objectnum,2) = objecttime; % destroy time
            end            
        end
        END_time = sscanf(str,'MSG %d'); % End Trial time = DestroyAll time
        for i=1:numel(trialobjects{trial})
            if lifetimes{trial}(i,2)==0
                lifetimes{trial}(i,2) = END_time;
            end
        end
%     elseif ~isempty(findstr(str,'Destroyed All Objects')) || ~isempty(findstr(str,'END TRIAL')) % Set any undestroyed objects to have their destroy time here
%         DA_time = sscanf(str,'MSG %d'); % DestroyAll time
%         for i=1:numel(trialobjects{trial})
%             if lifetimes{trial}(i,2)==0
%                 lifetimes{trial}(i,2) = DA_time;
%             end
%         end                
    end
end

                