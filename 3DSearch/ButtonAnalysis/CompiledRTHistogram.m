function [average, stdev, compiledRTVector] = CompiledRTHistogram (subjectNumber, inputFiles, initialTypeStr,binSize)
 
% CompiledRTHistogram
%
% Creates a histogram of the reaction times for ALL the mazes.
%
% [average, stdev, compiledRTVector] = CompiledRTHistogram
%    (subjectNumber, inputFiles, initialTypeStr,binSize)
% INPUT: subjectNumber is the subject's ID number in the data files to be
% loaded
% INPUT: inputFiles is a vector of the session numbers in the data files to
% be loaded
% INPUT: initialTypeStr is a string indicating which times should be used as
% anchors for reaction time
%   Input 'objects' for OBJECT APPEARANCE
%   Input 'saccades1' for FIRST SACCADE
%   Input 'saccades2' for LAST SACCADE
% INPUT: binSize is a scalar indicating how wide each bin should be in the
% histogram.
% OUTPUTS: average and stdev are the mean and standard deviation of the
% reaction times found
% OUTPUT: compiledRTVector is a vector of all the reaction times found.
% NOTE THAT ALL VALUES ARE IN UNITS OF MILLISECONDS
%
% Created 8/02/10 by Ansh Johri.
% Updated 8/03/10 by Ansh Johri - added user input for either saccades or 
% objects. Provides options for both FIRST saccade to object and SECOND
% saccade. Graph titles reflect type.
% Updated 2/15/11 by DJ - changed labels on histogram
% Updated 7/29/11 by DJ - binSize input

initialTypeNum = 0;

if nargin<4
    binSize = 50;
end

switch lower(initialTypeStr)        %% Determine the type of initial event. Display error if choice is unavailable.
    case {'objects'}
        initialTypeNum = 0;
    case {'saccades1'}
        initialTypeNum = 11;
    case {'saccades2'}
        initialTypeNum = 12;
    otherwise             
        error('Error.. please select "objects," "saccades1," or "saccades2." Defaulting to objects...')
end



compiledRTVector = [];
 
 for i = inputFiles  % Goes through each of the files
     fileName = sprintf('3DS-%d-%d',subjectNumber, i);  % Create file name
     load (fileName);  % Load file
     [a b curReactionTimes] = CreateRTHistogram(x, initialTypeStr);  % Get the reaction time vector for the particular file
     compiledRTVector = [compiledRTVector curReactionTimes];  % Add this vector into a vector comprised of ALL the reaction times
 end
 
average = mean(compiledRTVector);  % Calculate average reaction time
stdev = std(compiledRTVector);  % Calculate standard deviation
 
bin_centers = (binSize/2):binSize:max(compiledRTVector); % 'binSize' ms-wide bins from 0 to max RT
N = hist(compiledRTVector,bin_centers);  % Make histogram of this compiled vector, along with labels and title
bar(bin_centers,N,'hist');
xlabel('Reaction Time')

if initialTypeNum == 0  % Writes the title of the graph accordingly.
    titleLabel = 'Reaction Times of ALL Button Presses (Object Appearance)';
elseif initialTypeNum == 11
    titleLabel = 'Reaction Times of ALL Button Presses (First Saccade)';
elseif initialTypeNum == 12
    titleLabel = 'Reaction Times of ALL Button Presses (Last Saccade)';
end
    
meanLabel = sprintf('Mean: %g',average);
stdLabel = sprintf('STD: %g', stdev);
title({titleLabel;meanLabel; stdLabel})

