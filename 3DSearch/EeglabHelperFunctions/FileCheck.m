function filenameOut = FileCheck(filenameIn)

% Just a little gadget to check the desired eeglab filename.
%
% filenameOut = FileCheck(filenameIn)
%
% INPUT:
% -filenameIn is a string of the filename you want to check.
% OUTPUT:
% -filenameOut is a string of the ok'd filename in the current path.
%
% Created 2/17/11 by DJ.


fprintf('---Running FileCheck on %s...\n',filenameIn);

% Is it in the path?
if exist(filenameIn,'file')
    filenameOut = filenameIn;
    disp('File found!')
    return;


% If not...
else
    disp('File not found!')
    % Does it end in .set?
    setSpots = strfind(filenameIn,'.set');
    if isempty(setSpots)                % Filename lacks .set extension
        disp('Filename does not contain .set extention!')
        % Add .set extension
        newFilename = strcat(filenameIn,'.set');
        fprintf('Checking %s...\n',newFilename);
        
        if exist(newFilename,'file')    % Guess exists in current path
            disp('File found!')
            accept = input(sprintf('Suggested filename: %s.  Accept? (Y/N) >>',newFilename),'s'); % Ask user
            if strcmp(accept,'Y')       % User accepts guess
                filenameOut = newFilename;
                return;
            end
        else                            % Guess doesn't exist in current path
            disp('File not found!')
        end
        
    elseif length(setSpots)>1           % Filename contains too many .set extensions
        disp('Filename contains more than 1 .set extention!')
        % crop everything after first .set extension
        newFilename = filenameIn(1:(setSpots(1)+3)); % delete extra filename
        fprintf('Checking %s...\n',newFilename);
        
        if exist(newFilename,'file')    % Guess exists in current path
            disp('File found!')
            accept = input(sprintf('Suggested filename: %s.  Accept? (Y/N) >>',newFilename),'s'); % Ask user
            if strcmp(accept,'Y')       % User accepts guess
                filenameOut = newFilename;
                return;
            end
        else                            % Guess doesn't exist in current path
            disp('File not found!')    
        end
    end                                
        % At this point, the filename looks ok, but was rejected or not found in current path.
        
%     disp('Looking for other filenames for this subject...');    
%     % Are there other files for the same subject?
%     subject = sscanf(filenameIn,'%*s-%d-%*s');

    error('FileCheck could not predict proper filename.  Check for typos and add proper folders to path.')
    
end