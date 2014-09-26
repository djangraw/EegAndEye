function PlotPrePioBehavior(filename,pioThreshold,zEpoch,plottype)

% Plot the joystick input, path error, PIO amplitude, and pupil size just
% before the subject enters a PIO.
%
% PlotFlightSimBehavior(filename)
% 
% INPUTS:
% -filename is a string naming a file created using GetFlightSimBehavior.
% -pioThreshold is a scalar indicating the magnitude of the hilbert
% transform of the path error above which the subject is considered as
% entering a PIO. To lock to the end of the trial, set pioThreshold=inf.
% -zEpoch is a vector of z positions, relative to the point where the pio
% threshold is exceeded, at which the values should be recorded.
% -plottype is a string indicating the desired type of plot. 
%   + 'line': plot each individual trial as a line.
%   + 'density': plot a gray 2D density plot with a red line at the mean.
%
% Created 9/22/14 by DJ based on TestFlightSimFile.m.


if ~exist('zEpoch','var') || isempty(zEpoch)
    zEpoch = -2500:0;
end
if ~exist('plottype','var') || isempty(plottype)
    plottype = 'density';
end

% Load results, if they're already saved
if exist(filename,'file')
    fprintf('%s: Loading %s...\n',datestr(now,16),filename)
    load(filename);
    nTrials = numel(subjpos_cell);
else
    error('DJ:PlotFlightSimBehavior:FileNotFound','File %s not found!',filename);
end

% Set up
if ~isempty(pup_cell{1})
    nPlots = 4;
else
    nPlots = 3;
end


%% Get epochs
[ctrl_epoch, error_epoch, hilb_epoch, pup_epoch] = deal(nan(nTrials,length(zEpoch)));
zDone = zeros(1,nTrials);
for iTrial = 1:nTrials
    
    subjpos = subjpos_cell{iTrial};
    isInFlight = isInFlight_cell{iTrial};
    idealPath = [subjpos(isInFlight,1),interp1(ringpos(:,1),ringpos(:,2),subjpos(isInFlight,1),'linear','extrap')];
    pathError = [zeros(sum(~isInFlight),1); -(idealPath(:,2)-subjpos(isInFlight,2))];

    hilbPathError = abs(hilbert(pathError));
    hilbPathError(~isInFlight)=0;
    
    iDone = find(hilbPathError>pioThreshold,1);    
    if isempty(iDone)
        iDone = length(hilbPathError);
    end
    zDone(iTrial) = subjpos(iDone,1);
    zEpoch_this = zDone(iTrial) + zEpoch;
    
    % Fill epoch structs
    ctrl_epoch(iTrial,:) = interp1(ctrlPos_cell{iTrial}(:,1),abs(ctrl_cell{iTrial}),zEpoch_this);
    error_epoch(iTrial,:) = interp1(subjpos_cell{iTrial}(:,1),abs(pathError),zEpoch_this);    
    hilb_epoch(iTrial,:) = interp1(subjpos_cell{iTrial}(:,1),hilbPathError,zEpoch_this);    
    isInEpoch = ~isnan(pupPos_cell{iTrial});
    pup_epoch(iTrial,:) = interp1(pupPos_cell{iTrial}(isInEpoch),pup_cell{iTrial}(isInEpoch),zEpoch_this);
    
end

% subtract subject mean from all pupil size msmts
mean_pup = nanmean(cat(1,pup_cell{:}),1);
pup_epoch = pup_epoch - mean_pup;

zlimits = [min(zEpoch) max(zEpoch)];


%% Plot

n = 2^4; % resolution of plot

% --- Set up plot
fprintf('%s: Plotting Results...\n',datestr(now,16))
clf; 
subplot(nPlots,1,1); hold on;
ylimits = [0 2];
switch plottype
    case 'line'
        % plot subject controls
        plot(zEpoch,ctrl_epoch,'b');
    otherwise
        xdata = repmat(zEpoch,nTrials,1);
        ydata = ctrl_epoch;
        data = [xdata(:),ydata(:)];
        [~,density,X,Y] = kde2d(data,n,[zlimits(1) ylimits(1)],[zlimits(2) ylimits(2)]);
        % Plot
        cla; %hold on; % don't call 'hold on' until we get axis limits
        imagesc(X(1,:),Y(:,1),density);
        plot(zEpoch,nanmean(ctrl_epoch,1),'r','linewidth',2);

end

% Annotate plot
xlim(zlimits);
ylim(ylimits);
xlabel('Z position (m)');
ylabel('Joystick Input Amplitude (A.U.)');
title(show_symbols(filename))
% legend('Joystick Input','Ring Positions','Zero','Drift Change','Drift Amplitude','Location','EastOutside');





% --- Plot subject path
subplot(nPlots,1,2); hold on;

ylimits = [0 100];
switch plottype
    case 'line'
        % plot subject path error
        plot(zEpoch,error_epoch,'b');
    otherwise
        xdata = repmat(zEpoch,nTrials,1);
        ydata = error_epoch;
        data = [xdata(:),ydata(:)];
        [~,density,X,Y] = kde2d(data,n,[zlimits(1) ylimits(1)],[zlimits(2) ylimits(2)]);       
        
        % Plot
        imagesc(X(1,:),Y(:,1),density);
        plot(zEpoch,nanmean(error_epoch,1),'r','linewidth',2);

end


% Annotate plot
xlim(zlimits);
ylim(ylimits);
xlabel('Z position (m)');
ylabel('Elevation Error (m)');
% legend('Subject Path: takeoff','Subject Path: flight','Rings','Boundaries','Ground','Drift Change','Location','EastOutside');



% --- Set up PIO plot
subplot(nPlots,1,3); hold on;

ylimits = [0 100];
switch plottype
    case 'line'
        % plot PIO amplitude
        plot(zEpoch,hilb_epoch,'b');
    otherwise
        xdata = repmat(zEpoch,nTrials,1);
        ydata = hilb_epoch;
        data = [xdata(:),ydata(:)];
        [~,density,X,Y] = kde2d(data,n,[zlimits(1) ylimits(1)],[zlimits(2) ylimits(2)]);
        % Plot
        cla; %hold on; % don't call 'hold on' until we get axis limits
        imagesc(X(1,:),Y(:,1),density);
        plot(zEpoch,nanmean(hilb_epoch,1),'r','linewidth',2);
end


% Annotate plot
xlim(zlimits);
ylim(ylimits);
xlabel('Z position (m)');
ylabel('PIO Amplitude');
% ylabel('Hilbert Transform of error (m)');
% legend('hilb(error=0): takeoff','hilb(error): flight','Rings','Boundary','Drift Change','Location','EastOutside');


% --- Set up pupil size plot
if ~isempty(pupPos_cell{1})
    
    subplot(nPlots,1,4); hold on;
    
    ylimits = [-1000 1000];
    % Plot pupil size
    plot(zEpoch,pup_epoch,'b');
    switch plottype
        case 'line'
            % Plot pupil size
            plot(zEpoch,pup_epoch,'b');
        otherwise
            xdata = repmat(zEpoch,nTrials,1);
            ydata = pup_epoch;
            data = [xdata(:),ydata(:)];
            [~,density,X,Y] = kde2d(data,n,[zlimits(1) ylimits(1)],[zlimits(2) ylimits(2)]);
            % Plot
            cla; %hold on; % don't call 'hold on' until we get axis limits
            imagesc(X(1,:),Y(:,1),density);
            plot(zEpoch,nanmean(pup_epoch,1),'r','linewidth',2);
    end
    
    
    % Annotate plot
    xlim(zlimits);
    ylim(ylimits);
    xlabel('Z position (m)');
    ylabel('Pupil Size (A.U.)');
    % title(show_symbols(new_filename))
%     legend('Pupil Size','Location','EastOutside');
end


% Standardize plot widths
plotpos = zeros(nPlots,4);
for i=1:nPlots
    plotpos(i,:) = get(subplot(nPlots,1,i),'Position');
end
minwidth = min(plotpos(:,3));
plotpos(:,3) = minwidth;
for i=1:nPlots    
    set(subplot(nPlots,1,i),'Position',plotpos(i,:));    
end
colormap gray

fprintf('%s: Done!\n',datestr(now,16))
