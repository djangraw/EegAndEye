%% Load data
% EEG = pop_loadset('filename','3DS-22-all-filtered-noduds.set','filepath','/Users/dave/Documents/Data/3DSearch/2013-03-27-Pilot-S22/');
% load('3DS-22-2-eyepos');

%% Get data
EEGsamples = SmoothData(EEG.data(chans,(EEG.srate*12.45):(EEG.srate*14.95)),10);
EEGsamples(1,:) = EEGsamples(1,:) + 100;
EEGsamples2 = interp1((1:length(EEGsamples))/250,EEGsamples',(1:3000)/1000);

EYEsamples = SmoothData(eyepos(4850+(1:2500),:)',10)';

%% Set up movie
LIMITS = [275 730];

screen_res = [1024 768];
figure(111); clf;
% Make Eye Plot
h.EyePlot = axes('Units','normalized','Position',[0.13 0.4 0.775 0.55],'ydir','reverse','xtick',[],'ytick',[],'fontsize',16); % set position
title('Gaze');
hold on;
rectangle('Position',[0 0 screen_res],'linewidth',2);
plot([LIMITS(1) LIMITS(1)],[0 screen_res(2)],'g--','linewidth',2);
plot([LIMITS(2) LIMITS(2)],[0 screen_res(2)],'r--','linewidth',2);
h.Dot = plot(0,0,'k.','MarkerSize',50);
axis(h.EyePlot,[-200 screen_res(1)+200 -200 screen_res(2)+200]);
% Make EEG plot
h.EegPlot = axes('Units','normalized','Position',[0.13 0.1 0.775 0.23],'Ytick',[0 100],'Yticklabel',{'R','L'},'fontsize',16); % set position
title('EEG');
ylabel('Voltage')
xlabel('Time (ms)');
hold on;
plot([0 0; 2500 2500],[100 0; 100 0]);
plot(EEGsamples2,'linewidth',2);
ylims = get(gca,'ylim');
h.Line = plot([1 1],ylims,'k-','linewidth',2);
axis(h.EegPlot,[0 2500 -50 180 ]);


%% Start Movie
%     filename = uiputfile; % query user for filename
%     if isequal(filename,0)
%         return;
%     end

    fprintf(['Saving Movie to ' filename '...'])
    aviobj = avifile(filename,'fps',10); % slowed down 10x    




%% Play movie
newState = 0;
highlights = [];
for i=1:10:length(EYEsamples)
    oldState = newState;
    % Detect new state
    if EYEsamples(i,1)<LIMITS(1) 
        color = 'g';        
        newState = 1;       
    elseif EYEsamples(i,1)>LIMITS(2)
        color = 'r';
        newState = 2;
    elseif i>2050 && i<2300
        color = 'b';
        newState = 3;
    else
        color = 'k';
        newState = 0;
    end

    % move dot and line
    set(h.Dot,'xdata',EYEsamples(i,1),'ydata',EYEsamples(i,2),'color',color);
    set(h.Line,'xdata',[i i]);
    % Plot rectangle
    if newState>0 
        if newState==oldState
            oldrect = get(highlights{end},'xdata');
            set(highlights{end},'xdata',[oldrect(1:2);i;i]);
        else % add new highlight
            axes(h.EegPlot);
            highlights = {highlights patch([i;i;i;i],[ylims ylims([2 1])],color,'facealpha',0.5)};
        end
    end    
    pause(0.01);
    % Add to movie
    F = getframe(gcf);
    aviobj = addframe(aviobj,F);
end

%% Finish movie
    fprintf('\nClosing Movie...\n');
    aviobj = close(aviobj);
    disp('Done!')