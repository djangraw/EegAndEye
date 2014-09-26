function [EEG,components] = RemoveEogComponents(prefix,subject,sessions,eegSuffix,offset_ms)

% Removes HEOG and Blink components from a dataset with eye tracking data.
%
% EEG = RemoveEogComponents(prefix,subject,sessions,eegSuffix,offset_ms)
%
% INPUTS:
% -prefix is a string indicating the experiment code.
% -subject is a scalar indicating the subject number.
% -sessions is a vector indicating the session numbers.
% -eegSuffix is a string indicating the suffix on the EEG filename you want
% to use for both finding the components and removing the data. Example:
% '-filtered' (EEG file is <prefix>-<subject><eegSuffix>.set).
% -offset_ms is a scalar indicating the offset (in ms) of the reported 
% event times from the true times (as determined by matching up eye
% position with EOG).
%
% OUTPUTS:
% -EEG is an eeglab struct loaded from the specified file, with the HEOG
% and blink components removed.
%
% Created 6/5/13 by DJ.
% Updated 3/10/14 by DJ - added xMin option
% Updated 3/11/14 by DJ - added check to see if markers are already added

% Declare defaults
if nargin<5 || isempty(offset_ms)
    offset_ms = 0;
end
% select xMin automatically based on experiment type
if strcmp(prefix,'sq')
    xMin = 400;
else
    xMin = 200;
end
doPlot = true;

% Load data
disp('Loading...');
eegFilename = sprintf('%s-%d%s.set',prefix,subject,eegSuffix);
EEG = pop_loadset(eegFilename);
y = loadBehaviorData(subject,sessions,prefix);

% Check if events are already added
types = unique({EEG.event.type});

disp('Finding HEOG component...')
if any(strcmp('FSL',types))
    disp('Fixation markers already added...')
    EEG2 = EEG;
else
    EEG2 = AddHeogFixationMarkers(EEG,y,xMin);
end
comp_heog = GetHeogComponent(EEG2,offset_ms);

disp('Finding blink component...')
if any(strcmp('BS',types))
    disp('Blink markers already added...')
    EEG2 = EEG;
else
    EEG2 = AddBlinkMarkers(EEG,y);
end
comp_blink = GetBlinkComponent(EEG2,offset_ms);

disp('Removing components...')
components = [comp_heog, comp_blink];
EEG.data = SubtractOutComponents(EEG.data,components);

if doPlot
    disp('Plotting components...')
    figure;
    subplot(1,2,1);
    topoplot(double(comp_heog),EEG.chanlocs);
    title('HEOG component')
    subplot(1,2,2);
    topoplot(double(comp_blink),EEG.chanlocs);
    title('Blink component')
    MakeFigureTitle(EEG.setname);
end

disp('Done!')