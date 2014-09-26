function squares = MakeSquares(squareSize,dotSize,saveImages)

% Creates a set of Red/Green dot stimuli for the 3DSearch task.  Each
% stimulus is a black square with several pairs of colored rectangles in
% it.  Each colored rectangle is either red/green or green/red.
%
% squares = MakeSquares(squareSize,dotSize)
%
% INPUTS:
% - squareSize is the number of pixels on each side of the desired image.
% - dotSize is the number of each pixels on each side of the 'dot' (the
% red/green pair of recgangles.
% - saveImages is a binary value indicating whether the resulting images
% should be saved in the current directory (1) or plotted with imagesc (0).
%
% OUTPUTS:
% - squares is a 2^nDots element cell vector.  Each cell contains a 
% squareSize x squareSize x 3 matrix indicating an RGB image.
%
% Created 8/29/11 by DJ.

if nargin<3
    saveImages = false; % plot by default
end

% Set up
nDots = 4;
nPermutations = 2^nDots;
% isRG = zeros(nPermutations,nDots); % is dot red-on-left, green-on-right?

% Make RG and GR dots to paste into squares later
RG = zeros(dotSize,dotSize,3);
RG(:,1:dotSize/2,1) = 1; % red half
RG(:,dotSize/2+1:end,2) = 1; % green half
GR = zeros(dotSize,dotSize,3);
GR(:,1:dotSize/2,2) = 1; % green half
GR(:,dotSize/2+1:end,1) = 1; % red half

% Initialize variables
squares = cell(1,nPermutations);
isRG = zeros(1,nDots);

% Make squares
for i=1:nPermutations
    % Convert to binary values
    permstring = dec2bin(i-1,nDots);
    for j=1:length(permstring)
        isRG(j) = (permstring(j)=='1');
    end
    
    
    thisSquare = zeros(squareSize,squareSize,3); % RGB
    % top left square
    if isRG(1)
        thisSquare(round(squareSize/4-dotSize/2)+(1:dotSize),round(squareSize/4-dotSize/2)+(1:dotSize),:) = RG;        
    else
        thisSquare(round(squareSize/4-dotSize/2)+(1:dotSize),round(squareSize/4-dotSize/2)+(1:dotSize),:) = GR;  
    end
    % bottom left square
    if isRG(2)
        thisSquare(round(squareSize*3/4-dotSize/2)+(1:dotSize),round(squareSize/4-dotSize/2)+(1:dotSize),:) = RG;        
    else
        thisSquare(round(squareSize*3/4-dotSize/2)+(1:dotSize),round(squareSize/4-dotSize/2)+(1:dotSize),:) = GR;  
    end
    % top right square
    if isRG(3)
        thisSquare(round(squareSize/4-dotSize/2)+(1:dotSize),round(squareSize*3/4-dotSize/2)+(1:dotSize),:) = RG;        
    else
        thisSquare(round(squareSize/4-dotSize/2)+(1:dotSize),round(squareSize*3/4-dotSize/2)+(1:dotSize),:) = GR;  
    end
    % bottom right square
    if isRG(4)
        thisSquare(round(squareSize*3/4-dotSize/2)+(1:dotSize),round(squareSize*3/4-dotSize/2)+(1:dotSize),:) = RG;        
    else
        thisSquare(round(squareSize*3/4-dotSize/2)+(1:dotSize),round(squareSize*3/4-dotSize/2)+(1:dotSize),:) = GR;  
    end
    squares{i} = thisSquare;    
end

if saveImages
    for i=1:numel(squares)
        imwrite(squares{i},sprintf('square_%d.tiff',i),'tiff');
    end
    fprintf('Saved %d images to current directory\n',numel(squares));
else % plot them
    disp('press any key to advance to next stimulus...')
    for i=1:numel(squares)
        % Convert to binary values
        permstring = dec2bin(i-1,nDots);
        for j=1:length(permstring)
            isRG(j) = (permstring(j)=='1');
        end
        % plot
        imagesc(squares{i});
        title(['isRG = ' num2str(isRG)]);
        pause;
    end
    disp('Done!');
end