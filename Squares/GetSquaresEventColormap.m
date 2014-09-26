function Cmap = GetSquaresEventColormap(eventtypes)

% Get a colormap for consistent plotting of events.
%
% Cmap = GetSquaresEventColormap(legendnames)
%
% INPUTS:
% -eventtypes is a cell array of strings indicating the events in need of
% colors.
%
% OUTPUTS:
% - Cmap is an Nx3 matrix in which row i contains the RGB values for event
% type eventtypes{i}.
%
% Created 3/12/14 by DJ.
% Updated 3/18/14 by DJ - added SqNum indicators
% Updated 3/26/14 by DJ - added Square and Circle event types
% Updated 4/14/14 by DJ - changed D+ and T+ colors
% Updated 7/31/14 by DJ - added multiplier for ramp events

Cmap = zeros(length(eventtypes),3);
for j = 1:length(eventtypes)
    % unless we have reason to fade, make it bright.
    multiplier = 1; % factor by which to multiply this colormap
    % Find event type
    if ~isempty(strfind(eventtypes{j},'SqNum'))
        switch eventtypes{j}(end)
            case '1'
                Cmap(j,:) = [1 0 0]; %R
            case '2'
                Cmap(j,:) = [1 .5 0]; %O
            case '3'
                Cmap(j,:) = [1 1 0]; %Y
            case '4'
                Cmap(j,:) = [0 .75 0]; %G
            case '5'
                Cmap(j,:) = [0 0 1]; %B
            otherwise
        end
    elseif ~isempty(strfind(eventtypes{j},'Square'))
        Cmap(j,:) = [1 1 1]*.75; % LGy
    elseif ~isempty(strfind(eventtypes{j},'Circle'))
        Cmap(j,:) = [1 1 1]*.25; % DGy 

    else    
        
        if ~isempty(strfind(eventtypes{j},'-RampUp'))
            multiplier = 0.5;
            eventtypes{j} = eventtypes{j}(1:strfind(eventtypes{j},'-RampUp')-1);
        elseif ~isempty(strfind(eventtypes{j},'-RampDown'))
            multiplier = 0.5;
            eventtypes{j} = eventtypes{j}(1:strfind(eventtypes{j},'-RampDown')-1);
        elseif ~isempty(strfind(eventtypes{j},'-Peak'))
            multiplier = 0.5;
            eventtypes{j} = eventtypes{j}(1:strfind(eventtypes{j},'-Peak')-1);
        elseif ~isempty(strfind(eventtypes{j},'-Valley'))
            multiplier = 0.5;
            eventtypes{j} = eventtypes{j}(1:strfind(eventtypes{j},'-Valley')-1);
        else
            multiplier = 1;
        end
    
        switch eventtypes{j}(2:end-2) % exclude leading a/p and trailing 2/3
            case 'D_{0/'
                Cmap(j,:) = [0 0 1]; %B
            case 'D_{1/'
                Cmap(j,:) = [0 .75 .5]; %BG
            case {'D_{2/' 'D_{-/'}
                Cmap(j,:) = [0 1 1]; %C
            case 'D^*_{-/' 
                Cmap(j,:) = [0 .75 0]; %G
            case 'D_{+/'                        
                Cmap(j,:) = [0 .23 .63]; % Dark Blue                    
            case 'T_{0/'
                Cmap(j,:) = [1 0 0]; %R
            case {'T_{1/' 'T^*_{1/'}
                Cmap(j,:) = [.5 0 .75]; %V
            case {'T^*_{2/'}
                Cmap(j,:) = [1 0 1]; % M
            case 'T_{+/'
                Cmap(j,:) = [.63 0 .63]; % Dark Magenta
            case 'T_{-/'
                Cmap(j,:) = [1 1 1]*.25; % DGy 
            otherwise
                Cmap(j,:) = ([j j j]-1)/length(eventtypes); % Gy
        end
    end
    % Apply the multiplier to "fade" the event type
    Cmap(j,:) = Cmap(j,:)*multiplier;
end