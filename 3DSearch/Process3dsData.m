function [R, iObjects_eeg_pt] = Process3dsData(subject,sessions,offset,nSensitivity,cvmode)

% [R, iObjects_eeg_pt] = Process3dsData(subject,sessions,offset,nSensitivity,cvmode)
% 
% Created 4/2/13 by DJ as Process3dsData_script.
% Updated 5/13/13 by DJ - made into function.
% Updated 5/16/13 by DJ - added offset and nSensitivity inputs


if nargin<5 || isempty(cvmode)
    cvmode = '10fold';
end
%% Run ICA
% EEG = pop_runica(EEG, 'icatype','runica','dataset',1,'options',{'extended' 1});
% EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:end-4) '-ica.set'],'filepath',EEG.filepath);

%% Load ICA dataset
% subject = 22;
EEG = pop_loadset('filename',sprintf('3DS-%d-all-filtered-noduds-ica.set',subject));

%% Find eye components
% pop_topoplot(EEG,0, [1:20] ,EEG.setname,0 ,0,'electrodes','off'); % plot scalp maps
% pop_eegplot( EEG, 0, 1, 1); % plot component activations
% figure; pop_spectopo(EEG, 0, [0  EEG.pnts], 'EEG' , 'freq', [10], 'plotchan', 0, 'percent', 20, 'icacomps', [1:20], 'nicamaps', 5, 'freqrange',[2 25],'electrodes','off'); % plot spectra

%% Remove eye components
if subject==22
    BadICs = [1 3 12];
elseif subject==23
    BadICs = 1:4;
elseif subject==24
    BadICs = [1 2 4];
elseif subject==25
    BadICs = 1:3;
elseif subject==26
    BadICs = 2:4;
elseif subject==27
    BadICs = [3 4 7];
elseif subject==28
    BadICs = [2 3 5];
elseif subject==29
    BadICs = [2 3 5];
elseif subject==30
    BadICs = [1 4 7];
else
    % plot figures to help user decide on bad ICs
    pop_topoplot(EEG,0, [1:20] ,EEG.setname,0 ,0,'electrodes','off'); % plot scalp maps
    pop_eegplot( EEG, 0, 1, 1); % plot component activations
    figure; pop_spectopo(EEG, 0, [0  EEG.pnts], 'EEG' , 'freq', [10], 'plotchan', 0, 'percent', 20, 'icacomps', [1:20], 'nicamaps', 5, 'freqrange',[2 25],'electrodes','off'); % plot spectra
    % ask user for bad ICs
    BadICs = str2num(input('Input Bad Electrodes? >>','s'));
end
% BadICs = [1 2 4];
EEG = pop_subcomp( EEG, BadICs, 0);
BadIcString = sprintf('%d-',BadICs);
BadIcString(end) = [];
EEG.setname= sprintf('%s ICs_%s_Removed',EEG.setname, BadIcString);
EEG = eeg_checkset( EEG );

%% 20Hz LPF
% EEG = pop_iirfilt( EEG, 0, 20, [], [0]);
% EEG.setname= [EEG.setname, ' 20HzLPF'];

%% Epoch data
epochwin = [-1 2] + offset/1000; % in seconds
baselinewin = [0 100] + offset; % in ms
% baselinewin = [-100 0] + offset; % in ms

Numbers = GetNumbers;
sac2targ = num2str(Numbers.SACCADE_TO + Numbers.TARGET);
sac2dist = num2str(Numbers.SACCADE_TO + Numbers.DISTRACTOR);

EEG = pop_epoch( EEG, {sac2targ sac2dist}, epochwin, 'newname', [EEG.setname ' epochs'], 'epochinfo', 'yes');
EEG = pop_rmbase( EEG, baselinewin);% Remove baseline
EEG = eeg_checkset( EEG );

%% Enforce voltage threshold
vThresh = 50;
% EEG = EnforceVoltageThreshold(EEG,vThresh);
EEG = EnforceVoltageThreshold(EEG,vThresh,{sac2targ sac2dist},{},epochwin,[]);
EEG = pop_rejepoch(EEG,EEG.etc.rejectepoch,0);
EEG.setname = sprintf('%s %duVThresh',EEG.setname, vThresh);

%% Separate targets & distractors
clear ALLEEG
ALLEEG(1) = pop_epoch( EEG, { sac2dist  }, epochwin, 'newname', sprintf('%s DistractorEpochs',EEG.setname), 'epochinfo', 'yes');
ALLEEG(2) = pop_epoch( EEG, { sac2targ  }, epochwin, 'newname', sprintf('%s TargetEpochs',EEG.setname), 'epochinfo', 'yes');
% save(sprintf('TEMP_ALLEEG%s_old.mat',EEG.subject),'ALLEEG');
save(sprintf('ALLEEG%s_eyeposcorrected.mat',EEG.subject),'ALLEEG');

%% Run classifier on results
% cvmode = '10fold';
useica = true;

% get session #s
if nargin<2 || isempty(sessions)
    if subject==22
        sessions = 2:14;
    elseif subject==23
        sessions = [3 6:17];
    elseif subject>=24 && subject<=29
        sessions = 1:15;
    elseif subject==30
        sessions = [1:10 12:15];
    end
end
    
[iObjects_tag_pt, iObjects_eeg_pt,R] = RunResultsThroughTAG(subject,sessions,ALLEEG,cvmode,useica,offset,nSensitivity);
