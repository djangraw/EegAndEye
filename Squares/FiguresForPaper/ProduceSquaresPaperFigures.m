% ProduceSquaresPaperFigures
%
% Created 9/6/13 by DJ.
% Updated 3/19/14 by DJ - tResponse in cells

foo0 = load('sq-19-GLMresults-Theta_v2-LREvents-SqNum-Type-v2pt4CvsD0');
chanlocs = foo0.EEG.chanlocs;
if iscell(foo0.tResponse)
    tResponse = foo0.tResponse{foo0.iLevel};
else
    tResponse = foo0.tResponse;
end
% tBinCtr = 25:50:475;
tBinCtr = [125 275 425];
% tBinCtr = 250;
tBinWidth = 50;
clim = [-2 2];

%% Set plot params
% Set plot parameters
fontname = 'Helvetica'; % Futura, Helvetica, Gotham, Georgia are all good
fontsize = 10;
linewidth = 2;
markersize = 20;

%% Get RFs - all to 1
baseevent = 'C';
events = {'D0','D1','I','Inc'};
legendstr = {'Second Target vs.\n Early Distractor','Second Target vs.\n Late Distractor','Second Target vs.\n First Target','Second Target vs.\n Last Distractor'};

% baseevent = 'D0';
% events = {'D1','I','C','Inc'};
% legendstr = {'Late Distractor vs.\nEarly Distractor','First Target vs.\nEarly Distractor','Second Target vs.\nEarly Distractor','Last Distractor vs.\nEarly Distractor'};


RF = zeros(numel(chanlocs),numel(tBinCtr),numel(events));
isokevent = true(1,numel(events));
for i=1:numel(events)
    
    filename = sprintf('GroupGlm_sq_LREvents-SqNum-NewType-v2pt3%svs%s_2013-04-29.mat',events{i},baseevent);
    if exist(filename,'file')
        foo = load(filename);
        RF(:,:,i) = GetScalpMaps(foo.group_RF{1},tResponse,tBinCtr,tBinWidth);     
    else
        filename = sprintf('GroupGlm_sq_LREvents-SqNum-NewType-v2pt3%svs%s_2013-04-29.mat',baseevent,events{i});
        if exist(filename,'file')
            foo = load(filename);
            RF(:,:,i) = -GetScalpMaps(foo.group_RF{1},tResponse,tBinCtr,tBinWidth);     
        else
            fprintf('%s vs. %s comparison does not exist! Filling with NaNs.\n',events{i},baseevent);
            RF(:,:,i) = NaN;
            isokevent(i) = false;
        end
    end
end

%% Get RFs - pairwise
% eventpairs = {'IvsD0','CvsI','CvsD0','CvsD1','CvsInc','D1vsD0','IncvsD1'};
% legendstr = {'T_0_/_2 - D_0_/_2','T_1_/_2^* - T_0_/_2','T_1_/_2^* - D_0_/_2',...
%     'T_1_/_2^* - D_1_/_2','T_1_/_2^* - D_1_/_2^*','D_1_/_2 - D_0_/_2',...
%     'D_1_/_2^* - D_1_/_2'};

% eventpairs = {'D1vsD0','IvsD0','CvsD0'};
% legendstr = {'D_1_/_2 - D_0_/_2','T_0_/_2 - D_0_/_2','T_1_/_2^* - D_0_/_2'};
% fLegendstr = {'fD_1_/_2 - fD_0_/_2','fT_0_/_2 - fD_0_/_2','fT_1_/_2^* - fD_0_/_2'};
% sLegendstr = {'sD_1_/_2 - sD_0_/_2','sT_0_/_2 - sD_0_/_2','sT_1_/_2^* - sD_0_/_2'};

eventpairs = {'IvsD0', 'CvsD0', 'CvsI'};
legendstr = {'T_0_/_2 - D_0_/_2', 'T_1_/_2^* - D_0_/_2', 'T_1_/_2^* - T_0_/_2'};
fLegendstr = {'fT_0_/_2 - fD_0_/_2', 'fT_1_/_2^* - fD_0_/_2', 'fT_1_/_2^* - fT_0_/_2'};
sLegendstr = {'sT_0_/_2 - sD_0_/_2', 'sT_1_/_2^* - sD_0_/_2', 'sT_1_/_2^* - sT_0_/_2'};

% eventpairs = {'IvsD0', 'CvsD1', 'CvsI'};
% legendstr = {'T_0_/_2 - D_0_/_2', 'T_1_/_2^* - D_1_/_2', 'T_1_/_2^* - T_0_/_2'};
% fLegendstr = {'fT_0_/_2 - fD_0_/_2', 'fT_1_/_2^* - fD_1_/_2', 'fT_1_/_2^* - fT_0_/_2'};
% sLegendstr = {'sT_0_/_2 - sD_0_/_2', 'sT_1_/_2^* - sD_1_/_2', 'sT_1_/_2^* - sT_0_/_2'};

% eventpairs = {'S2vsS1','S3vsS1','S4vsS1','S5vsS1','CirvsS1'};
% legendstr = {'Sq2 - Sq1','Sq3 - Sq1','Sq4 - Sq1','Sq5 - Sq1','Cir - Sq1'};
% fLegendstr = {'fSq2 - fSq1','fSq3 - fSq1','fSq4 - fSq1','fSq5 - fSq1','Cir - fSq1'};
% sLegendstr = {'sSq2 - sSq1','sSq3 - sSq1','sSq4 - sSq1','sSq5 - sSq1','Cir - sSq1'};

% eventpairs = {'S1vsS5','S2vsS5','S3vsS5','S4vsS5'};
% legendstr = {'Sq1 - Sq5','Sq2 - Sq5','Sq3 - Sq5','Sq4 - Sq5'};
% fLegendstr = {'fSq1 - fSq5','fSq2 - fSq5','fSq3 - fSq5','fSq4 - fSq5'};
% sLegendstr = {'sSq1 - sSq5','sSq2 - sSq5','sSq3 - sSq5','sSq4 - sSq5'};


% legendstr = {'First Target vs.\nEarly Distractor','Second Target vs.\nEarly Distractor',...
%     'Second Target vs.\nLate Distractor','Second Target vs.\nFirst Target',...
%     'Last Distractor vs.\nLate Distractor','Second Target vs.\nLast Distractor'};

% bandnames = {'', 'Theta_v2-' 'Lower1Alpha_v2-' 'Lower2Alpha_v2-' 'UpperAlpha_v2-','Delta_v2-' 'Beta_v2-' 'LowerGamma_v2-'};
% bandnames = {'Theta_v2-' 'Lower1Alpha_v2-' 'Lower2Alpha_v2-' 'UpperAlpha_v2-','Beta_v2-' 'LowerGamma_v2-'};
% bandlimits = [6, 8, 10, 12, 14, 28, 48];
bandnames = {''};
[sRF,sZ,fRF,fZ] = deal(zeros(numel(chanlocs),numel(tResponse),numel(eventpairs),numel(bandnames)));
% bandname = 'Theta_v2-';

for iBand = 1:numel(bandnames)
    bandname = bandnames{iBand};

    for i=1:numel(eventpairs)

        % Saccade
    %     filename = sprintf('GroupGlm_sq_%sLREvents-SqNum-NewType-v2pt3%s_2013-04-29.mat',bandname,eventpairs{i});
        sFilename = sprintf('GroupGlm_sq_%sLREvents-SqNum-Type-v2pt4%s_*.mat',bandname,eventpairs{i});    
    %     sFilename = sprintf('GroupGlm_sq_%sLREvents-Type-SqNum-v2pt4%s_2013-09-20.mat',bandname,eventpairs{i});    
        sFileFound = dir(sFilename);
        if ~isempty(sFileFound)        
            sFoo = load(sFileFound.name);                
            sRF(:,:,i,iBand) = sFoo.group_RF{1};     
            sZ(:,:,i,iBand) = sFoo.group_z{1};
        else
            fprintf('sq, %s: %s comparison does not exist! Filling with NaNs.\n',bandnames{iBand},eventpairs{i});
            sRF(:,:,i,iBand) = NaN;
            sZ(:,:,i,iBand) = NaN;
        end

        % Fixation
        fFilename = sprintf('GroupGlm_sf_%sSqNum-Type-v2pt4%s_*.mat',bandname,eventpairs{i});
%         fFilename = sprintf('GroupGlm_sf_%sSqNum-Type-v2pt4%s_2013-09-12.mat',bandname,eventpairs{i});
    %     fFilename = sprintf('GroupGlm_sf_%sSqNum-Type-v2pt4%s_2013-09-13.mat',bandname,eventpairs{i});
    %     fFilename = sprintf('GroupGlm_sf_%sType-SqNum-v2pt4%s_2013-09-20.mat',bandname,eventpairs{i});
        fFileFound = dir(fFilename);    
        if ~isempty(fFileFound)         
            fFoo = load(fFileFound.name);               
            fRF(:,:,i,iBand) = fFoo.group_RF{1};
            fZ(:,:,i,iBand) = fFoo.group_z{1};             
        else
            fprintf('sf, %s: %s comparison does not exist! Filling with NaNs.\n',bandnames{iBand},eventpairs{i});
            fRF(:,:,i,iBand) = NaN;
            fZ(:,:,i,iBand) = NaN;
        end

    end
end
%% Plot RFs
thingsToPlot = {'fRF','fZ','sRF','sZ'};
iBand = 1;

for i=1:numel(thingsToPlot)
    option = thingsToPlot{i}(2:end);
    eval(sprintf('thingToPlot = %s(:,:,:,iBand);',thingsToPlot{i}));

    RF = GetScalpMaps(thingToPlot,tResponse,tBinCtr,tBinWidth);
    isokevent = ~isnan(RF(1,1,:));
    RF_ok = RF(:,:,isokevent);
    legendstr_ok = legendstr(isokevent);

    switch option    
        case 'RF'
            if isempty(bandname)
                clim = [-2 2];
            else
                clim = [-.2 .2];
            end
            cthresh = 0;
        case 'Z';
            clim = [-5 5];
            cthresh = 1.96;
    end

    figure(840+i); clf;
    % set(gcf,'Position',[1 480 2560 1000]);
    PlotScalpMaps(RF_ok,chanlocs,clim,'','',cthresh);
end

%% make labels
for i=1:4
    figure(840+i);
    option = thingsToPlot{i}(2:end);
    switch option    
        case 'RF'
            if isempty(bandname)
                clim = [-2 2];
            else
                clim = [-.2 .2];
            end
            cthresh = 0;
        case 'Z';
            clim = [-5 5];
            cthresh = 1.96;
    end
    
    for j=1:size(RF_ok,2)
        subplot(size(RF_ok,3),size(RF_ok,2),size(RF_ok,2)*(size(RF_ok,3)-1)+j);
        xlabel(sprintf('%d',tBinCtr(j)),'visible','on','fontname',fontname,'fontsize',fontsize);
    end
    for j = 1:size(RF_ok,3)
        subplot(size(RF_ok,3),size(RF_ok,2),(j-1)*size(RF_ok,2)+1);
        ylabel(sprintf(legendstr_ok{j}),'visible','on','fontname',fontname,'fontsize',fontsize);
    end

    subplot(size(RF_ok,3),size(RF_ok,2),size(RF_ok,2)*size(RF_ok,3));
    lowerright = get(gca,'Position');
    axispos = [lowerright(1)+lowerright(3)/2, lowerright(2:4)];
    axes('Position',axispos,'CLim',clim,'visible','off');%,'fontname',fontname,'fontsize',fontsize);
    colormap(ThresholdedCmap(clim,cthresh));
    colorbar
end

%% Plot RF timecourses

chanGrid = {'F3','FZ','F4'; 'C3','CZ','C4'; 'P3','PZ','P4'};
thingToPlot = fRF;
% thingToPlot = cat(3,fRF,sRF);
option = 'RF';
legendstr_ok = legendstr;
% legendstr_ok = [fLegendstr, sLegendstr];

figure(948); clf;
% set(gcf,'Position',[1 480 2560 1000]);
PlotResponseFnsGrid(thingToPlot,legendstr_ok,tResponse,chanlocs,chanGrid,Cmap);
if strcmp(option,'Z');
    for i=1:numel(chanGrid)
        subplot(size(chanGrid,1),size(chanGrid,2),i);
        hold on;
        plot(get(gca,'xlim'),[cthresh, cthresh],'k--','linewidth',2);
        plot(get(gca,'xlim'),[-cthresh, -cthresh],'k--','linewidth',2);
    end
end

%% Plot SVD decomp
thingToPlot = cat(3,fRF,sRF);
% legendstr_ok = legendstr;
legendstr_ok = [fLegendstr, sLegendstr];
tRange = [0 500];
nComps = 4;
sigma = 0;

figure(950);clf;
thingToPlot_smooth = SmoothData(thingToPlot,sigma);
[weights, eigenvalues] = ApplySvdToGlmResults(thingToPlot_smooth,tResponse,tRange,chanlocs,legendstr_ok,nComps);

%% Plot occipital electrode average
elecsToAvg = { 'PO9'    'PO7'    'PO3'    'PO1'    'POZ'    'PO2'    'PO4'    'PO6'    'PO8'    'PO10'...
    'O1'    'OZ'    'O2'};
fakelabel = struct('labels','Occipital Avg');
% elecsToAvg = { 'C5'    'C3'    'C1'    'CZ'    'C2'    'C4'    'C6'   ...
%     'CP5'    'CP3'    'CP1'    'CPZ'    'CP2'    'CP4'    'CP6'   'P5'    'P3'    'P1'    'PZ'...
%     'P2'    'P4'    'P6'};
% fakelabel = struct('labels','Centro-Parietal Avg');
iElecs = find(ismember({chanlocs.labels},elecsToAvg));
sigma = 5; % in samples

thingToPlot = cat(3,mean(sZ(iElecs,:,:),1));%,mean(sRF(iElecs,:,:),1));
thingToPlot_smooth = SmoothData(thingToPlot,sigma);

legendstr_ok = [fLegendstr];%,sLegendstr];
figure(962); clf;
% set(gcf,'Position',[1 480 2560 1000]);
PlotResponseFns(thingToPlot_smooth,legendstr_ok,tResponse,fakelabel,1);

