function EEG = eogartifacts(EEG,FIX_EXPAND,OUTSIDE_WINDOW)

% Given an unepoched EEGLAB data file with saccade events 
% included, find the maximum power difference component between L-R
% saccades and R-L saccades, and project it out of that saccade data.
% 
% EEG = eogartifacts(EEG,FIX_EXPAND,OUTSIDE_WINDOW)
%
% INPUTS:
% - EEG is an unepoched EEGLAB data file with saccade events
% included. (Currently the 3DSearch codes for these events are used.)
% - FIX_EXPAND is a scalar indicating how far eog artifacts extend 
% beyond the boundaries defined by eyelink events, or two scalars
% indicating the time relative to the start and end fixation events that
% should be considered part of the saccade (e.g. [-50, 100]).
% - OUTSIDE_WINDOW is a scalar indicating how much data on either side of
% each fixation should be averaged during mean subtraction.  if this = inf,
% the component is removed from the entire dataset.
%
% OUTPUTS:
% - EEG is the same dataset but with the saccade component removed from the
% data inside saccade events, and ' - NoEogArtifacts' added to EEG.setname.
%
% Created 10/28/11 based on blinkartifacts.
% Updated 11/28/11 to work with text events

if nargin<2
    FIX_EXPAND = 0; % time, in ms, beyond the saccade start/end events surrounding a saccade, that you want to be considered in the saccade.
end
if nargin<3
    OUTSIDE_WINDOW = 25; % time, in ms, on either side of the saccade that we want to average to do mean subtraction
end

if length(FIX_EXPAND)==1
    FIX_EXPAND = [-FIX_EXPAND FIX_EXPAND];
end

disp('---Removing HEOG Artifacts:')
disp('Removing HEOG Component...')
% Get time info
% extract event info
eventTypes = str2double({EEG.event(:).type});
eventLatencies = [EEG.event(:).latency]' * 1000/EEG.srate;

% Detect file type
if strcmp(EEG.setname(1:2),'sq')
    Constants = GetSquaresConstants;
    if mean(isnan(eventTypes))<0.5 % Mostly numeric codes
        left_start_times = eventLatencies(ismember(eventTypes,Constants.FIXSTART_BASE+Constants.LEFTSIDE));
        left_end_times = eventLatencies(ismember(eventTypes,Constants.FIXEND_BASE+Constants.LEFTSIDE));
        right_start_times = eventLatencies(ismember(eventTypes,Constants.FIXSTART_BASE+Constants.RIGHTSIDE));
        right_end_times = eventLatencies(ismember(eventTypes,Constants.FIXEND_BASE+Constants.RIGHTSIDE));            
    else % Mostly string codes
        eventTypes = {EEG.event(:).type};
        left_start_times = eventLatencies(strmatch(Constants.EVENTNAMES{Constants.FIXSTART_BASE+Constants.LEFTSIDE},eventTypes));
        left_end_times = eventLatencies(strmatch(Constants.EVENTNAMES{Constants.FIXEND_BASE+Constants.LEFTSIDE},eventTypes));
        right_start_times = eventLatencies(strmatch(Constants.EVENTNAMES{Constants.FIXSTART_BASE+Constants.RIGHTSIDE},eventTypes));
        right_end_times = eventLatencies(strmatch(Constants.EVENTNAMES{Constants.FIXEND_BASE+Constants.RIGHTSIDE},eventTypes));            
    end
else
    error('Only squares files are currently supported!')
end
leftRanges = [left_start_times left_end_times];
rightRanges = [right_start_times right_end_times];

% expand the saccade ranges by 25ms in either direction
leftRanges(:,1) = leftRanges(:,1)+FIX_EXPAND(1);
leftRanges(:,2) = leftRanges(:,2)+FIX_EXPAND(2);
rightRanges(:,1) = rightRanges(:,1)+FIX_EXPAND(1);
rightRanges(:,2) = rightRanges(:,2)+FIX_EXPAND(2);

% Determine whether each point is within a saccade
times = (1:size(EEG.data,2))/EEG.srate*1000; % times in ms
isInLFix = false(1,size(EEG.data,2));
isInRFix = false(1,size(EEG.data,2));
for i=1:size(leftRanges,1)
    isInLFix(times>=leftRanges(i,1) & times<=leftRanges(i,2)) = true;        
end
for i=1:size(rightRanges,1)
    isInRFix(times>=rightRanges(i,1) & times<=rightRanges(i,2)) = true;
end
isInFix = isInLFix | isInRFix;

% Concatenate saccade data
leftData = EEG.data(:,isInLFix);
rightData = EEG.data(:,isInRFix);
saccadeData = EEG.data(:,isInFix);

% Find max power component
A = difference(leftData, rightData);

% Get vectors for eye subtraction
V = (A'*A)\A'; % Pseudoinverse of A

% Remove component from saccade areas
cleanFixData = (eye(size(saccadeData,1))-A*V)*saccadeData; % x_new = x - (v^T*x)/(v^T*v)

% Add clean data back to EEG struct
if isinf(OUTSIDE_WINDOW) % special case outlined in header
    % apply to entire dataset
    EEG.data = (eye(size(saccadeData,1))-A*V)*EEG.data;
%     return;
else

    % Place clean data back in data matrix
    EEG.data(:,isInFix) = cleanFixData;

    % --- Perform mean subtraction 
    if OUTSIDE_WINDOW>0 % if there is any outside data for mean subtraction
        disp('Performing Mean Subtraction...');    
        for i=1:size(leftRanges,1)
            thisFix = (times>=leftRanges(i,1) & times<=leftRanges(i,2));
            if (leftRanges(i,2)-leftRanges(i,1)) < (2*OUTSIDE_WINDOW) % short saccade --> include all inside data
                thisInside = thisFix;
            else % long saccade --> only use data near edges
                thisInside = ( ( times>=leftRanges(i,1) & times<leftRanges(i,1)+OUTSIDE_WINDOW ) | ( times>leftRanges(i,2)-OUTSIDE_WINDOW & times<=leftRanges(i,2) ) );
            end
            thisOutside = ( ( times>=leftRanges(i,1)-OUTSIDE_WINDOW & times<leftRanges(i,1) ) | ( times>leftRanges(i,2) & times<=leftRanges(i,2)+OUTSIDE_WINDOW ) );
            meanDiff = mean(EEG.data(:,thisInside),2)-mean(EEG.data(:,thisOutside),2);
            EEG.data(:,thisFix) = EEG.data(:,thisFix) - repmat(meanDiff,1,sum(thisFix));    
        end
    end
    
end

% % --- Check for discontinuities
% discontinuity = zeros(nFixs,2);
% for i=1:nFixs
%     firstInside = find(times>=leftRanges(i,1),1,'first');    
%     lastInside = find(times<=leftRanges(i,2),1,'last');
%     discontinuity(i,1) = mean(abs(EEG.data(:,firstInside-1)-EEG.data(:,firstInside)));
%     discontinuity(i,2) = mean(abs(EEG.data(:,lastInside)-EEG.data(:,lastInside+1)));
% end
% figure;
% subplot(1,2,1);
% hist(discontinuity(:,1),100);
% xlabel('discontinuity before saccade (uV)')
% ylabel('number of saccades')
% subplot(1,2,2);
% hist(discontinuity(:,2),100);
% xlabel('discontinuity after saccade (uV)')
% ylabel('number of saccades')
% fom = mean(discontinuity(:)); % figure of merit

% --- Update EEG struct
% Change setname to mark that we've done this
EEG.setname=[EEG.setname ' - NoEOGArtifacts'];

disp('Plotting...')
% --- Plot results
figure;
% Plot component
subplot(1,2,1)
topoplot(A,EEG.chanlocs,'electrodes','on');
colorbar
title('Max Difference EOG Component')
% Find frontal electrodes
plotLabels = {'Fp1' 'Fpz' 'Fp2'};
chanLabels = {EEG.chanlocs(:).labels}; % the names of all the channels
iChansToPlot = find(ismember(chanLabels,plotLabels));
if isempty(iChansToPlot)
    plotLabels = {'FP1','FPZ','FP2'};
    iChansToPlot = find(ismember(chanLabels,plotLabels));
end

% Plot data before and after
firstFixLength = find(times<leftRanges(1,2)-leftRanges(1,1),1,'last'); % number of samples in first saccade
subplot(2,2,2)
plot(times(1:firstFixLength),saccadeData(iChansToPlot,1:firstFixLength) + repmat(50*(1:length(iChansToPlot))',1,firstFixLength));
title('First left-side fixation: Before')
xlabel('time (ms)')
ylabel([sprintf('voltage on electrodes\n') sprintf('%s, ',plotLabels{:}) 'in uV'])
ylimits = get(gca,'ylim');
subplot(2,2,4)
plot(times(1:firstFixLength),cleanFixData(iChansToPlot,1:firstFixLength) + repmat(50*(1:length(iChansToPlot))',1,firstFixLength));
title('First left-side saccade: After')
xlabel('time (ms)')
ylabel([sprintf('voltage on electrodes\n') sprintf('%s, ',plotLabels{:}) 'in uV'])
set(gca,'ylim',ylimits)
MakeFigureTitle(sprintf('RemoveEogArtifacts ([%g %g], %g)',FIX_EXPAND(1),FIX_EXPAND(2),OUTSIDE_WINDOW));

disp('Done!')


% --- Find the component representing the difference between the means in 
% the given data.
% This function was taken from eyesubtract.m.
function Difvec=difference(EEG_1,EEG_2)
    % Find component
    Difvec = mean(EEG_1,2)-mean(EEG_2,2);
    % Normalize
    Difvec   = Difvec./norm(Difvec);

return