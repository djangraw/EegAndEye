function EEG = blinkartifacts(EEG,BLINK_EXPAND,OUTSIDE_WINDOW)

% Given an unepoched EEGLAB data file with blink and saccade events 
% included, find the maximum power component within all data inside blinks
% and project it out of that blink data.
% 
% EEG = blinkartifacts(EEG,BLINK_EXPAND,OUTSIDE_WINDOW)
%
% INPUTS:
% - EEG is an unepoched EEGLAB data file with blink and saccade events
% included. (Currently the 3DSearch codes for these events are used.)
% - BLINK_EXPAND is a scalar indicating how far blink artifacts extend 
% beyond the boundaries defined by eyelink events, or two scalars
% indicating the time relative to the start and end blink events that
% should be considered part of the blink (e.g. [-50, 100]).
% - OUTSIDE_WINDOW is a scalar indicating how much data on either side of
% each blink should be averaged during mean subtraction.  if this = inf,
% the component is removed from the entire dataset.
%
% OUTPUTS:
% - EEG is the same dataset but with the blink component removed from the
% data inside blink events, and ' - NoBlinkArtifacts' added to EEG.setname.
%
% Created 8/11/11 by DJ (named RemoveBlinkArtifacts).
% Updated 8/16/11 by DJ - changed to blinkartifacts, made
% RemoveBlinkArtifacts a wrapper, added BLINK_EXPAND 2-element option.
% Updated 10/26/11 by DJ - added all-caps labels for Sensorium-2011.
% Updated 11/23/11 by DJ - added hard-coded kluge for sq-6-all dataset.

if nargin<2
    BLINK_EXPAND = 0; % time, in ms, beyond the saccade start/end events surrounding a blink, that you want to be considered in the blink.
end
if nargin<3
    OUTSIDE_WINDOW = 25; % time, in ms, on either side of the blink that we want to average to do mean subtraction
end

if length(BLINK_EXPAND)==1
    BLINK_EXPAND = [-BLINK_EXPAND BLINK_EXPAND];
end

disp('---Removing Blink Artifacts:')
disp('Removing Blink Component...')
% Get time info
blinkRanges = BlinkRange(EEG); % blink start and end times in ms
times = (1:size(EEG.data,2))/EEG.srate*1000; % times in ms

% expand the blink ranges by 25ms in either direction
blinkRanges(:,1) = blinkRanges(:,1)+BLINK_EXPAND(1);
blinkRanges(:,2) = blinkRanges(:,2)+BLINK_EXPAND(2);

% Determine whether each point is within a blink
nBlinks = size(blinkRanges,1);
isInBlink = false(1,size(EEG.data,2));
for i=1:nBlinks
    isInBlink(times>=blinkRanges(i,1) & times<=blinkRanges(i,2)) = true;
end
if strncmp(EEG.setname,'sq-6-all',8)
    isInBlink(times>2170e3) = 0; % KLUGE to get rid of crazy filter waves at end of last session!
end

% Concatenate blink data
blinkData = EEG.data(:,isInBlink);
blinkData(blinkData>1000) = 0;

% Find max power component
A = Maximumpower(blinkData);

% Get vectors for eye subtraction
V = (A'*A)\A'; % Pseudoinverse of A

% Remove component from blink areas
cleanBlinkData = (eye(size(blinkData,1))-A*V)*blinkData; % x_new = x - (v^T*x)/(v^T*v)

% Add clean data back to EEG struct
if isinf(OUTSIDE_WINDOW) % special case outlined in header
    % apply to entire dataset
    EEG.data = (eye(size(blinkData,1))-A*V)*EEG.data;
    return;
else

    % Place clean data back in data matrix
    EEG.data(:,isInBlink) = cleanBlinkData;

    % --- Perform mean subtraction 
    if OUTSIDE_WINDOW>0 % ifC there is any outside data for mean subtraction
        disp('Performing Mean Subtraction...');    
        for i=1:nBlinks
            thisBlink = (times>=blinkRanges(i,1) & times<=blinkRanges(i,2));
            if (blinkRanges(i,2)-blinkRanges(i,1)) < (2*OUTSIDE_WINDOW) % short blink --> include all inside data
                thisInside = thisBlink;
            else % long blink --> only use data near edges
                thisInside = ( ( times>=blinkRanges(i,1) & times<blinkRanges(i,1)+OUTSIDE_WINDOW ) | ( times>blinkRanges(i,2)-OUTSIDE_WINDOW & times<=blinkRanges(i,2) ) );
            end
            thisOutside = ( ( times>=blinkRanges(i,1)-OUTSIDE_WINDOW & times<blinkRanges(i,1) ) | ( times>blinkRanges(i,2) & times<=blinkRanges(i,2)+OUTSIDE_WINDOW ) );
            meanDiff = mean(EEG.data(:,thisInside),2)-mean(EEG.data(:,thisOutside),2);
            EEG.data(:,thisBlink) = EEG.data(:,thisBlink) - repmat(meanDiff,1,sum(thisBlink));    
        end
    end
    
end

% % --- Check for discontinuities
% discontinuity = zeros(nBlinks,2);
% for i=1:nBlinks
%     firstInside = find(times>=blinkRanges(i,1),1,'first');    
%     lastInside = find(times<=blinkRanges(i,2),1,'last');
%     discontinuity(i,1) = mean(abs(EEG.data(:,firstInside-1)-EEG.data(:,firstInside)));
%     discontinuity(i,2) = mean(abs(EEG.data(:,lastInside)-EEG.data(:,lastInside+1)));
% end
% figure;
% subplot(1,2,1);
% hist(discontinuity(:,1),100);
% xlabel('discontinuity before blink (uV)')
% ylabel('number of blinks')
% subplot(1,2,2);
% hist(discontinuity(:,2),100);
% xlabel('discontinuity after blink (uV)')
% ylabel('number of blinks')
% fom = mean(discontinuity(:)); % figure of merit

% --- Update EEG struct
% Change setname to mark that we've done this
EEG.setname=[EEG.setname ' - NoBlinkArtifacts'];

disp('Plotting...')
% --- Plot results
figure;
% Plot component
subplot(1,2,1)
topoplot(A,EEG.chanlocs,'electrodes','on');
colorbar
title('Max Power Blink Component')
% Find frontal electrodes
plotLabels = {'Fp1' 'Fpz' 'Fp2'};
chanLabels = {EEG.chanlocs(:).labels}; % the names of all the channels
iChansToPlot = find(ismember(chanLabels,plotLabels));
if isempty(iChansToPlot)
    plotLabels = {'FP1','FPZ','FP2'};
    iChansToPlot = find(ismember(chanLabels,plotLabels));
end

% Plot data before and after
firstBlinkLength = find(times<blinkRanges(1,2)-blinkRanges(1,1),1,'last'); % number of samples in first blink
subplot(2,2,2)
plot(times(1:firstBlinkLength),blinkData(iChansToPlot,1:firstBlinkLength) + repmat(50*(1:length(iChansToPlot))',1,firstBlinkLength));
title('First blink: Before')
xlabel('time (ms)')
ylabel([sprintf('voltage on electrodes\n') sprintf('%s, ',plotLabels{:}) 'in uV'])
ylimits = get(gca,'ylim');
subplot(2,2,4)
plot(times(1:firstBlinkLength),cleanBlinkData(iChansToPlot,1:firstBlinkLength) + repmat(50*(1:length(iChansToPlot))',1,firstBlinkLength));
title('First blink: After')
xlabel('time (ms)')
ylabel([sprintf('voltage on electrodes\n') sprintf('%s, ',plotLabels{:}) 'in uV'])
set(gca,'ylim',ylimits)
MakeFigureTitle(sprintf('RemoveBlinkArtifacts ([%g %g], %g)',BLINK_EXPAND(1),BLINK_EXPAND(2),OUTSIDE_WINDOW));

disp('Done!')


% --- Find the component with maximum power (sum squared activity about the
% mean) in the given data.
% This function was taken from eyesubtract.m.
function Max_eigvec=Maximumpower(EEGdata)

[Channel,~,Epoch]=size(EEGdata);

% Subtract the mean from each channel... should we do this?
for e=1:Epoch
    for i=1:Channel
        EEGdata(i,:,e)=EEGdata(i,:,e)-mean(EEGdata(i,:,e));
    end;
end;

% Maximum Power method
a=EEGdata(:,:);

pow_a=a*a';

[vb,tmp] = eig(pow_a); % the second output is necessary here - it changes the first output results.
Max_eigvec=vb(:,end)./norm(vb(:,end));

return