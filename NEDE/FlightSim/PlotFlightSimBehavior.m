function PlotFlightSimBehavior(filename)

% PlotFlightSimBehavior(filename)
% 
% INPUTS:
%
% Created 9/22/14 by DJ based on TestFlightSimFile.m.


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
    nPlots = 5;
else
    nPlots = 4;
end


%% Plot
% --- Set up plot
fprintf('%s: Plotting Results...\n',datestr(now,16))
clf; 
subplot(nPlots,1,1); hold on;
ylim([-4 4]);
% Plot subject path
plot(ctrlPos_cell{1}(:,1),SmoothData(ctrl_cell{1},5,'half'));
% Plot ideal path (through ring centers)
plot(ringpos(:,1),zeros(size(ringpos,1),1),'r.-');
% plot ground
plot(zlimits,[0 0],'k--');

% plot rest of trials
for iTrial=1:nTrials    
    plot(ctrlPos_cell{iTrial}(:,1),SmoothData(ctrl_cell{iTrial},5,'half'));
end

% Annotate plot
xlim(zlimits);
xlabel('Z position (m)');
ylabel('Joystick Input (A.U.)');
title(show_symbols(filename))
legend('Joystick Input','Ring Positions','Zero','Drift Change','Drift Amplitude','Location','EastOutside');





% --- Plot subject path
subplot(nPlots,1,2); hold on;
subjpos = subjpos_cell{1};
isInFlight = isInFlight_cell{1};
plot(subjpos(~isInFlight,1),subjpos(~isInFlight,2),'b--');
plot(subjpos(isInFlight,1),subjpos(isInFlight,2));
% Plot ideal path (through ring centers)
plot(ringpos(:,1),ringpos(:,2),'r.-');
% Plot top edge
plot(ringpos(:,1),ringpos(:,2)+ringsize(:,2)/2,'r:');

% plot bottom edge (here to fix legend)
plot(ringpos(:,1),ringpos(:,2)-ringsize(:,2)/2 ,'r:');
% Plot rings
for i=1:size(ringpos,1)
    plot([ringpos(i,1), ringpos(i,1)], [ringpos(i,2)-ringsize(i,2)/2, ringpos(i,2)+ringsize(i,2)/2],'r-');        
end
% plot rest of trials
for iTrial=1:nTrials
    subjpos = subjpos_cell{iTrial};
    isInFlight = isInFlight_cell{iTrial};
    plot(subjpos(~isInFlight,1),subjpos(~isInFlight,2),'b--');
    plot(subjpos(isInFlight,1),subjpos(isInFlight,2));
end

% Annotate plot
xlim(zlimits);
xlabel('Z position (m)');
ylabel('Elevation (m)');
title(show_symbols(filename))
legend('Subject Path: takeoff','Subject Path: flight','Rings','Boundaries','Ground','Drift Change','Location','EastOutside');



% --- Set up PIO plot
subplot(nPlots,1,3); hold on;
subjpos = subjpos_cell{1};
isInFlight = isInFlight_cell{1};
idealPath = [subjpos(isInFlight,1),interp1(ringpos(:,1),ringpos(:,2),subjpos(isInFlight,1))];
pathError = [zeros(sum(~isInFlight),1); -(idealPath(:,2)-subjpos(isInFlight,2))];

hilbPathError = hilbert(pathError);

% plot(subjpos(:,1),abs(pathError));
plot(subjpos(~isInFlight,1),abs(hilbPathError(~isInFlight)),'b--');
plot(subjpos(isInFlight,1),abs(hilbPathError(isInFlight)));
plot(ringpos(:,1),zeros(size(ringpos,1),1),'r.-');
plot(ringpos(:,1),ringsize(:,2)/2,'r:');
    
% plot individual rings
for i=1:size(ringpos,1)
    plot([ringpos(i,1), ringpos(i,1)], [0 ringsize(i,2)/2],'r-');        
end

% plot rest of trials
for iTrial = 1:nTrials
    subjpos = subjpos_cell{iTrial};
    isInFlight = isInFlight_cell{iTrial};
    idealPath = [subjpos(isInFlight,1),interp1(ringpos(:,1),ringpos(:,2),subjpos(isInFlight,1),'linear','extrap')];
    pathError = [zeros(sum(~isInFlight),1); -(idealPath(:,2)-subjpos(isInFlight,2))];

    hilbPathError = hilbert(pathError);

    % plot(subjpos(:,1),abs(pathError));
    plot(subjpos(~isInFlight,1),abs(hilbPathError(~isInFlight)),'b--');
    plot(subjpos(isInFlight,1),abs(hilbPathError(isInFlight)));
end

xlim(zlimits);
xlabel('Z position (m)');
ylabel('PIO Amplitude');
% ylabel('Hilbert Transform of error (m)');
legend('hilb(error=0): takeoff','hilb(error): flight','Rings','Boundary','Drift Change','Location','EastOutside');


% --- Set up Missed Ring histogram
subplot(nPlots,1,4); cla;
ringZ = ringpos(:,1);
zEnd = zeros(1,nTrials);
for iTrial=1:nTrials
    zEnd(iTrial) = subjpos_cell{iTrial}(end,1);
end
nMissed = hist(zEnd,ringZ); 
plot(ringZ,nMissed);

xlim(zlimits);
xlabel('Z position (m)');
ylabel('Rings Missed');
legend('Histogram of missed rings','Location','EastOutside');



% --- Set up pupil size plot
if ~isempty(pupPos_cell{1})
    
    subplot(nPlots,1,5); hold on;
    % Plot subject path
    pup = pup_cell{1};
    pupPos = pupPos_cell{1};
    plot(pupPos,pup-mean(pup));
    % Plot ideal path (through ring centers)
    plot(ringpos(:,1),zeros(size(ringpos,1),1),'r.-');
    % plot ground
    plot(get(gca,'xlim'),[0 0],'k--');

    for iTrial = 1:nTrials
        pup = pup_cell{iTrial};
        pupPos = pupPos_cell{iTrial};
        plot(pupPos,pup-mean(pup));
    end
    % Annotate plot
    xlim(zlimits);
    xlabel('Z position (m)');
    ylabel('Pupil Size (A.U.)');
    % title(show_symbols(new_filename))
    legend('Pupil Size','Ring Positions','Zero','Drift Change','Location','EastOutside');
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

fprintf('%s: Done!\n',datestr(now,16))
