% function RunAbmNedeData(subject,sessions,ascPrefix,eegPrefix,eegFilename,tpFilename)
%%
subject = 101;
sessions = 2:11;
ascPrefix = 'NEDE_0pt3';
eegPrefix = 'NEDE';
rawFilename = 'Raw_Data.csv';
tpFilename = 'TP_Data.csv';

% basedir = '/Users/dave/Documents/Data/NEDE/ABM_Pilots'; % Tesla
basedir = '/Users/jangrawdc/Documents/LiincData/NEDE/ABM_Pilots'; % Dave's NIH machine
set(0,'DefaultFigureColormap',feval('jet')); % override MATLAB2014b's default

%%
cd(basedir)
[tData,eegData,tEvents,eegEvents] = deal(cell(1,numel(sessions)));
for i=1:numel(sessions)
    fprintf('%s: Importing session %d/%d...\n',datestr(now,16),i,numel(sessions));
    %%    
    ascFilename = sprintf('%s-%d-%d',ascPrefix,subject,sessions(i));
    NEDE_ImportData(ascPrefix,subject,sessions(i));
    
    % Import EEG data
    eegFilename = sprintf('%s_%d_%d/%s',eegPrefix,subject,sessions(i),rawFilename);
    values = csvread(eegFilename,1,0);
    fid = fopen(eegFilename,'r');    
%     headers = textscan(fid,repmat('%s',1,size(values,2)),1,'Delimiter',',');
%     headers = [headers{3:13}];
    fclose(fid);
    
    tData{i} = values(:,1);
    eegData{i} = values(:,3:11); % timestamps, EKG, EEG.
    
    % Import Third-Party events
    eventsFilename = sprintf('%s_%d_%d/%s',eegPrefix,subject,sessions(i),tpFilename);
    values = csvread(eventsFilename,1,0);
    fid = fopen(eventsFilename,'r');    
%     headers = textscan(fid,repmat('%s',1,size(values,2)),1,'Delimiter',',');    
%     headers = [headers{:}];
    fclose(fid);
    
    tEvents{i} = values(:,1);
    eegEvents{i} = values(:,2);   
    
end
%% Add Sync Events
raw_fs = 256;
new_fs = raw_fs;
output_suffix = '-filtered';
bandpass_bounds = [1 100];
notch_bounds = [59 61];
electrodeloc_filename = 'dave_abm9chan.loc';

for i=1:numel(sessions)
    
    % Import data 
    fprintf('%s: Importing session %d/%d to EEGLAB...\n',datestr(now,16),i,numel(sessions));
    EEG = NEDE_FilterEegData(eegData{i}',ascPrefix,subject,sessions(i), raw_fs,new_fs, output_suffix, bandpass_bounds,notch_bounds); 
    
    % Add events and re-save
    times = (tEvents{i}-tData{i}(1))/1000;
    codes = eegEvents{i};
    events = [times, codes]; % the times (in s) and codes of each event
    assignin('base','events',events);
    EEG = pop_importevent( EEG, 'append','yes','event','events','fields',{'latency' 'type'},'timeunit',1,'optimalign','off');
    
    EEG = pop_saveset(EEG,'filename',EEG.filename);
end

%% Combine datasets

input_suffix = '-filtered';
output_suffix = '-all-filtered';

EEG = CombineEeglabSessions(subject,sessions,input_suffix,output_suffix,ascPrefix);

%% Add events
% EEG = pop_loadset(sprintf('%s-%s%s.set',ascPrefix,subject,output_suffix));
new_filename = sprintf('%s-%d-all-events.set',ascPrefix,subject);

clear y
for i=1:numel(sessions)
    load(sprintf('%s-%d-%d.mat',ascPrefix,subject,sessions(i)));
    y(i) = x;
end
[eyeTimes, eyeCodes] = NEDE_GetEeglabEvents(y,EEG);
EEG = NEDE_AddEeglabEvents(y,EEG,eyeTimes,eyeCodes);
EEG = pop_saveset(EEG,'filename',new_filename);


%% Remove HEOG and Blink Components
% Add artifact event markers
EEG = pop_loadset(new_filename);
EEG = pop_reref(EEG,[]); % NEW 11/7!!! RE-REFERENCE TO AVG.
EEG = NEDE_AddHeogEvents(EEG,y);
% Remove components
blinktypes = {'Blink Start','Blink End'};
heogtypes = {'FSL','FEL';'FSR','FER'};
offset_ms = 0;
blinkComp = GetBlinkComponent(EEG,offset_ms,blinktypes);
heogComp = GetHeogComponent(EEG,offset_ms,heogtypes);
% Subtract out components
EEG.data = SubtractOutComponents(EEG.data,[blinkComp,heogComp]);

%% Plot results
figure;
subplot(1,2,1);
topoplot(double(blinkComp),EEG.chanlocs);
title('Blink compoenent');
subplot(1,2,2);
topoplot(double(heogComp),EEG.chanlocs);
title('HEOG compoenent');
MakeFigureTitle(new_filename);

%% Run ICA and save
EEG = pop_runica(EEG, 'icatype','runica','dataset',1,'options',{'extended' 1});
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:end-4) '-reref-ica.set'],'filepath',EEG.filepath);

%% Check ICs
ica_filename = sprintf('%s-%d-all-events-reref-ica.set',ascPrefix,subject);
EEG = pop_loadset(ica_filename);
nComps = size(EEG.icaweights,1);
pop_topoplot(EEG,0, 1:nComps ,EEG.setname,0 ,0,'electrodes','on'); % plot scalp maps
pop_eegplot( EEG, 0, 1, 1); % plot component activations
figure; pop_spectopo(EEG, 0, [0  EEG.pnts], 'EEG' , 'freq', [10], 'plotchan', 0, 'percent', 20, 'icacomps', 1:nComps, 'nicamaps', 5, 'freqrange',[2 25],'electrodes','off'); % plot spectra

%% Remove bad ICs
BadICs = 8:9;
EEG = pop_subcomp( EEG, BadICs, 0);
BadIcString = sprintf('%d-',BadICs);
BadIcString(end) = [];
EEG.setname= sprintf('%s ICs_%s_Removed',EEG.setname, BadIcString);
EEG = eeg_checkset( EEG );
% EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);

%% Epoch data
EEGtarg = pop_epoch(EEG,{'Targ Saccade'},[-1 2]);
% EEGtarg = pop_rmbase(EEGtarg,[0 100]); % remove post-saccade baseline
EEGdist = pop_epoch(EEG,{'Dist Saccade'},[-1 2]);
% EEGdist = pop_rmbase(EEGdist,[0 100]); % remove post-saccade baseline

nComps = size(EEG.icaweights,1);
EEGtarg.icaact = nan(nComps, EEGtarg.pnts, EEGtarg.trials);
for i=1:EEGtarg.trials
    EEGtarg.icaact(:,:,i) = (EEG.icaweights*EEG.icasphere)*EEGtarg.data(EEG.icachansind,:,i);
end

EEGdist.icaact = nan(nComps, EEGdist.pnts, EEGdist.trials);
for i=1:EEGdist.trials
    EEGdist.icaact(:,:,i) = (EEG.icaweights*EEG.icasphere)*EEGdist.data(EEG.icachansind,:,i);
end

%% Train classifier
twl_ms = 100; % in ms
two_ms = 100:100:900; % in ms
two = round(interp1(EEGtarg.times, 1:EEGtarg.pnts, two_ms));% indices
twl = round(twl_ms/1000*EEGtarg.srate);
%% Set up data
% data = cat(3,EEGtarg.data,EEGdist.data);
data = cat(3,EEGtarg.icaact(1:7,:,:),EEGdist.icaact(1:7,:,:));
fmdata = cat(3,EEGtarg.data,EEGdist.data);
truth = [ones(1,EEGtarg.trials), zeros(1,EEGdist.trials)];


%% Run HDCA classifier
% [y_level2, w, v, fwdModel, y_level1] = TrainHybridHdcaClassifier(data,truth,twl,two,[],fmdata,'nocrossval');
[y_level2, w, v, fwdModel, y_level1] = RunHybridHdcaClassifier(data,truth,twl,two,'10fold',[],fmdata);
%% Run Old RSVP Classifier
useica = true;
[y_level2, w, v, fwdModel] = run_rsvp_classifier([EEGdist EEGtarg],twl,two,'20fold',true);
% [y_level2, w, v, fwdModel] = run_rsvp_classifier([EEGdist EEGtarg],twl,two,'loo',useica);
%% Plot classifier results

PlotHybridHdcaClassifier(fwdModel, v, EEG.chanlocs, two_ms+twl_ms/2);

%%
clf;
Az_level1 = zeros(1,length(two));
for i=1:length(two);
    Az_level1(i) = rocarea(y_level1(:,i),truth);
end
plot(two_ms+twl_ms/2,Az_level1);