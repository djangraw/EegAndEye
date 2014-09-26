function PlotFixationProperties(x,samples)

% Plots information about fixations over the time an object is visible.
%
% PlotFixationProperties(x)
%
% -
%
% Created 9/13/10 by DJ.
% Updated 9/14/10 by DJ - 


% SET UP
% Set up visible times info
visibletimes = GetVisibleTimes(x,'eyelink');
obj = visibletimes(:,1);
isTarget = strcmp('TargetObject',{x.objects(obj).tag});
appeartimes = visibletimes(:,2);
disappeartimes = visibletimes(:,3);        
object_limits = x.eyelink.object_limits;

% Set up plot
clf;
max_y1 = 500; % distance from object in pixels
max_y2 = 1000; % duration in ms
max_t = 1000; % time after object appearance to be viewed
iPlot = 0; % which subplot we're currently on

% MAIN CODE
fix_times = x.eyelink.fixation_times(:,2); % fixation END time
fix_positions = x.eyelink.fixation_positions;
fix_durations = x.eyelink.fixation_times(:,2)-x.eyelink.fixation_times(:,1);
sample_times = x.eyelink.record_time-1 + (1:length(samples));
         
% Find all fixations made during each appearance and find distance from obj
for j=1:numel(appeartimes)
    
    %---FIXATION STUFF---
    % Find the distance from each applicable fixation to the object's bounding box
    fix_subset = find(fix_times>appeartimes(j) & fix_times < disappeartimes(j));
    fix_dist = Inf(numel(fix_subset),1); % the distance to the bounding box (preallocate for speed)
    t_fix = fix_times(fix_subset); % times when these saccades were made
    for k=1:numel(fix_subset)
        % Get the time when this object was last seen
        iLastLimit = find(object_limits(:,2)==obj(j) & object_limits(:,1)<t_fix(k),1,'last');
        % Get object limits
        left = object_limits(iLastLimit,3);
        top = object_limits(iLastLimit,4);
        width = object_limits(iLastLimit,5);
        height = object_limits(iLastLimit,6);
        % Get fixation position
        px = fix_positions(fix_subset(k),1);
        py = fix_positions(fix_subset(k),2);
        % Get shortest distance from fixation position to object box
        if px<left
            dx = left-px;
        elseif px>left+width
            dx = px-(left+width);
        else 
            dx = 0;
        end
        if py<top
            dy = top-py;
        elseif py>top+height
            dy = py-(top+height);
        else
            dy = 0;
        end
        fix_dist(k) = sqrt(dx*dx + dy*dy);
    end
    
    % ---SAMPLE STUFF---
    sample_subset = find(sample_times>appeartimes(j) & sample_times < disappeartimes(j));
    sample_dist = Inf(numel(sample_subset),1); % the distance to the bounding box (preallocate for speed)
    t_sample = sample_times(sample_subset); % times when these saccades were made
    
    for l=1:numel(sample_subset)
        % Get the time when this object was last seen
        iLastLimit = find(object_limits(:,2)==obj(j) & object_limits(:,1)<t_sample(l),1,'last');
        % Get object limits
        left = object_limits(iLastLimit,3);
        top = object_limits(iLastLimit,4);
        width = object_limits(iLastLimit,5);
        height = object_limits(iLastLimit,6);
        % Get fixation position
        px = samples(sample_subset(l),1);
        py = samples(sample_subset(l),2);
        % Get shortest distance from fixation position to object box
        if px<left
            dx = left-px;
        elseif px>left+width
            dx = px-(left+width);
        else 
            dx = 0;
        end
        if py<top
            dy = top-py;
        elseif py>top+height
            dy = py-(top+height);
        else
            dy = 0;
        end
        sample_dist(l) = sqrt(dx*dx + dy*dy);
    end
    
    % ---OTHER INFO---
    % get the duration of each fixation
    dur = fix_durations(fix_subset); % the distance to the bounding box (preallocate for speed)
    % get the fraction of the object that was visible
    iLims = (object_limits(:,2)==obj(j) & object_limits(:,1)>appeartimes(j) & object_limits(:,1)<disappeartimes(j));
    fracvis_x = object_limits(iLims,1) - appeartimes(j);
    fracvis_y = object_limits(iLims,7) * max_y1;
    
    % ---PLOT---
    %Plot lines over time
    if ~isempty(t_fix)
        iPlot = iPlot+1;
        if iPlot>30 break; end
        subplot(5,6,iPlot);
        ax = plotyy(t_fix-appeartimes(j),fix_dist,t_fix-appeartimes(j),dur);
        hold on;
        plot(fracvis_x,fracvis_y,'r');
        plot(t_sample-appeartimes(j),sample_dist,'c');
        
        xlabel('time from object appearance (ms)')
        title(sprintf('%d: Obj %d, %s',appeartimes(j),obj(j),x.objects(obj(j)).tag));
        set(ax(1),'Ylim',[0 max_y1],'YTickMode','auto','XLim',[0 max_t]);
        set(ax(2),'Ylim',[0 max_y2],'YTickMode','auto','XLim',[0 max_t]);
        set(get(ax(1),'Ylabel'),'String','distance (pix)');
        set(get(ax(2),'Ylabel'),'String','fix duration (ms)');
    end
end
MakeLegend({'b','c','g','r'},{'fixpos','samplepos','fixdur','frac vis'});
MakeFigureTitle(sprintf('Fixation Stats for Subject %d, Session %d',x.subject,x.session));

% xlabel('time from object appearance (ms)')
% ylabel('distance from target (pixels)');
