function [average, stdev, singleRTVector] = CreateRTHistogram(mainStruct, initialTypeStr)

% Creates a histogram of the reaction times for the particular subject and session.
%
% [average, stdev, singleRTVector] = CreateRTHistogram(mainStruct, initialTypeNum)
%
% INPUTS: 
% -mainStruct is the main structure file containing the particular session.
% -initialTypeStr is a string indicating which times should be used as
% anchors for reaction time.
%   Input 'objects' for OBJECT APPEARANCE
%   Input 'saccades1' for FIRST SACCADE
%   Input 'saccades2' for LAST SACCADE
% OUTPUTS: 
% -average and stdev are the mean and standard deviation of the
% reaction times found
% -compiledRTVector is a vector of the reaction times found for this
% particular session.
% NOTE THAT ALL VALUES ARE IN UNITS OF MILLISECONDS
%
% Created 8/02/10 by Ansh Johri.
% Updated 8/03/10 by Ansh Johri - added user input for either saccades or 
% objects. Provides options for both FIRST saccade to object and SECOND saccade. 
% Updated 2/15/11 by DJ - now button_times can be a column vector

initialTypeNum = 0;


switch lower(initialTypeStr)        %% Determine the type of initial event. Display error if choice is unavailable.
    case {'objects'}
        initialTypeNum = 0;
    case {'saccades1'}
        initialTypeNum = 11;
    case {'saccades2'}
        initialTypeNum = 12;
    otherwise             
        error('Error.. please select "objects," "saccades1," or "saccades2."')
end


outputArray = CheckButtonPress(mainStruct, initialTypeStr);  % Initialize variables
singleRTVector = [];
samplingFactor = mainStruct.eeg.eventsamplerate/ 1000;

y = length(mainStruct.eeg.button_times);

for i = 1:y  % Loop for the each button press
    isTarget = outputArray(1,i);
    if isTarget == 1  || initialTypeNum == 11 || initialTypeNum == 12 % If button is pressed on a target object..
        timeTarget = outputArray(2,i);
        if (isnan(timeTarget))
        else
            timeButton = mainStruct.eeg.button_times(i)/samplingFactor;   % Get time that the button was pressed
            timeReaction = timeButton - timeTarget;  % Reaction time is the difference between these two times
            singleRTVector = [singleRTVector, timeReaction];  % Store reaction time in a vector
        end
    end
end

average = mean(singleRTVector);  % Calculate average reaction time
stdev = std(singleRTVector);  % Calculate standard deviation

hist(singleRTVector);  % Make histogram with labels and title
xlabel('Reaction Time');
title('Reaction Times of Button Presses');