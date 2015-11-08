% ImportAllSquaresScript
%
% Created 2/20/14 by DJ.
% Updated 2/19/15 by DJ - use GetSquaresSubjects

experiments = {'SquaresFix'};%{'Squares','SquaresFix','SquaresFix3'};
driftcorrection = true;
filetype = 'Sensorium-2013';
pixelthresh = 75;
bandpassbounds = [1 40];
notchbounds = [0 -1];
newfs = 100;
includeeog = false;
singlesuffix = '-40Hz-fs100';
combosuffix = '-all-40Hz-fs100';

for iExp = 1:numel(experiments)
    experiment = experiments{iExp};
    [subjects, basedir, folders] = GetSquaresSubjects(experiment);
    switch experiment
        case 'SquaresFix3'
            prefix = 'sf3';            
%             subjects = 1:12;
            eyeSessions_cell = {1:10, 2:11, 1:10, 1:10, 1:10, 1:10, 1:10, [1 2 5:12], 1:10, 1:10, 1:10, 1:10};
            eegSessions_cell = {2:11, 2:11, 0:9, 0:9, 0:9, 0:9, 0:9, [0 1 6:13], 0:9, 0:9, 0:9, 0:9};

%             basedir = '/Users/dave/Documents/Data/SquaresFix3';
%             folders = {'2013-05-15-Pilot-S1'    '2013-12-02-Pilot-S2'    '2013-12-04-Pilot-S3'    '2013-12-04-Pilot-S4',...
%                 '2013-12-05-Pilot-S5'    '2013-12-16-Pilot-S6'    '2013-12-19-Pilot-S7'    '2014-02-04-Pilot-S8',...
%                 '2014-02-06-Pilot-S9'    '2014-02-06-Pilot-S10'   '2014-02-12-Pilot-S11'   '2014-02-18-Pilot-S12'};

        case 'SquaresFix'
            prefix = 'sf';
%             subjects = [1:10 12:13];
            eyeSessions_cell = {1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, [1:5 7:11], 1:10};
            eegSessions_cell = {0:9, 0:9, 0:9, 0:9, 0:9, 0:9, 0:9, 0:9, [0:3 6:11], 0:9, [0:4 6:10], 0:9};

%             basedir = '/Users/dave/Documents/Data/SquaresFix';
%             folders = {'2013-03-05-Pilot-S1','2013-03-11-Pilot-S2','2013-03-18-Pilot-S3',...
%                 '2013-03-19-Pilot-S4','2013-04-01-Pilot-S5','2013-04-02-Pilot-S6',...
%                 '2013-04-09-Pilot-S7','2013-04-11-Pilot-S8','2013-04-11-Pilot-S9',...
%                 '2013-04-17-Pilot-S10','2014-02-17-Pilot-S12','2014-03-06-Pilot-S13'};
 
        case 'Squares'
            prefix = 'sq';
            
%             subjects = [9:11 13:15 17:19 20:27];
            eyeSessions_cell = {2:11, 1:10, 2:11, 1:10, 1:10, 1:10, 1:10, 2:11, [1 3:11], 1:10, 1:10, [1:6 8:11], 1:10, 1:10, 1:10,1:10,4:13};
            eegSessions_cell = {1:10, 0:9, 2:11, 0:9, 0:9, 0:9, 0:9, 1:10, [0 2:10], 0:9, [0 2:10], [0:5 7:10], 0:9, 0:9, 0:9,0:9,3:12};


%             folders = {'2012-03-29-Pilot-S9', '2012-04-09-Pilot-S10', '2012-04-18-Pilot-S11', ...
%                 '2012-05-04-Pilot-S13', '2012-05-29-Pilot-S14', '2012-06-12-Pilot-S15', ...
%                 '2012-06-15-Pilot-S17', '2012-06-21-Pilot-S18', '2012-06-22-Pilot-S19',...
%                 '2012-06-25-Pilot-S20', '2013-04-15-Pilot-S21', '2013-04-16-Pilot-S22', ...
%                 '2014-02-13-Pilot-S23', '2014-02-17-Pilot-S24', '2014-02-20-Pilot-S25', ...
%                 '2014-03-04-Pilot-S26', '2014-03-10-Pilot-S27'};

%             basedir = '/Users/dave/Documents/Data/Squares';

    end


    eventsrule = [experiment '-basic'];

    for iSubj = 1:numel(subjects)
        cd(basedir);
        cd(folders{iSubj});
        subject = subjects(iSubj);
        sessions = eyeSessions_cell{iSubj};
        eegsessions = eegSessions_cell{iSubj};

        if exist(sprintf('%s-%d%s.set',prefix,subject,combosuffix),'file')~=0
            fprintf('Skipping %s-%d import...\n',prefix,subject);
            continue;
        end
        
        
%         % Remove old files, if they're there
%         delete('*fs125*');
        
        % FILE CHECK!
        eeg_filename = sprintf('%s-%d.%.3d.bh',prefix,subject,eegsessions(1)); % .dat EEG file
        if exist([eeg_filename '.dat'],'file')~=0
            filetype = 'Sensorium-2013';
        else 
            filetype = 'Sensorium-2011';            
        end               
        
        for i=1:numel(sessions)
            % Run functions as in ImportData
            switch experiment
                case 'Squares'
                    import_squares_data(subject,sessions(i),eegsessions(i),filetype,pixelthresh,driftcorrection);
                case {'SquaresFix' 'SquaresFix3'}
                    import_squaresfix_data(subject,sessions(i),eegsessions(i),filetype,pixelthresh,driftcorrection,prefix);
            end
            ImportToEeglab_2014(subject,sessions(i),eegsessions(i),filetype,experiment,singlesuffix,bandpassbounds,notchbounds,newfs,includeeog);
            AddEeglabEvents(subject,sessions(i),singlesuffix,eventsrule);
        end
        %Combine sessions
        CombineEeglabSessions(subject,sessions,singlesuffix,combosuffix,experiment);
        saveBehaviorData(subject,sessions,prefix);
    end
end


%% Interpolate Bad Electrodes

combosuffix_interp = [combosuffix '-interp'];

outerElecs = {'NZ','F9','F10','FT9','FT10','T9','T10','TP9','TP10',...
    'P9','P10','PO9','PO10','I1','I2','IZ'};


for iExp = 1:numel(experiments)
    experiment = experiments{iExp};
    [subjects,basedir,folders] = GetSquaresSubjects(experiment);
    switch experiment
        case 'SquaresFix3'
            prefix = 'sf3';
%             subjects = 1:12;            
            badElecs_cell = {{'PO5','T9','T10','TP9','FP1','I2'},...
                {'T10'},...
                {'T9','T10'},...
                {'PO5','T9','T10','PO7'},...
                {'FC4'},... % T10 a little
                {'T10','TP9'},...
                {'TP9'},... % 'FP1','FPZ' early on
                {'T10'},...
                {'T9','T10','TP10'},... % they're not that bad
                {'PO5','T9','T10','T7'},...
                {'PO5','T9','T10','T8','TP9','CP4','TP8'},...
                {}};
            
%             basedir = '/Users/dave/Documents/Data/SquaresFix3';
%             folders = {'2013-05-15-Pilot-S1'    '2013-12-02-Pilot-S2'    '2013-12-04-Pilot-S3'    '2013-12-04-Pilot-S4',...
%                 '2013-12-05-Pilot-S5'    '2013-12-16-Pilot-S6'    '2013-12-19-Pilot-S7'    '2014-02-04-Pilot-S8',...
%                 '2014-02-06-Pilot-S9'    '2014-02-06-Pilot-S10'   '2014-02-12-Pilot-S11'   '2014-02-18-Pilot-S12'};

        case 'SquaresFix'
            prefix = 'sf';
%             subjects = [1:10 12:13];
            badElecs_cell = {{'PO5','FP1','T9','T10','I2'}, ...
                {'PO5','FP1','TP9','P9','I2'}, ...
                {'PO5','FP1','T9','T7','T10','TP7','PO9','I2'}, ...
                {'PO5','FP1','T9','T10','CPZ','I2'}, ...
                {'PO5','FP1','T9','T7','T10','TP9','I1','I2'}, ...
                {'PO5','FP1','T9','T10','CPZ','I2'}, ...
                {'PO5','FP1','T9','T10','CPZ','P9','P7','I2'}, ...
                {'PO5','FP1','T9','T10','CPZ','I2'}, ...
                {'PO5','FP1','T9','T10','CPZ','I2'}, ...
                {'PO5','FP1','T9','T10','CPZ','I2'}, ...
                {}, ... % CP6 is noisy after ~1820/2660 s
                {} }; %13
            
%             basedir = '/Users/dave/Documents/Data/SquaresFix';
%             folders = {'2013-03-05-Pilot-S1','2013-03-11-Pilot-S2','2013-03-18-Pilot-S3',...
%                 '2013-03-19-Pilot-S4','2013-04-01-Pilot-S5','2013-04-02-Pilot-S6',...
%                 '2013-04-09-Pilot-S7','2013-04-11-Pilot-S8','2013-04-11-Pilot-S9',...
%                 '2013-04-17-Pilot-S10','2014-02-17-Pilot-S12','2014-03-06-Pilot-S13'};
 
        case 'Squares'
            prefix = 'sq';            
%             subjects = [9:11 13:15 17:19 20:27];
            badElecs_cell = {{'PO5','T9','T10','TP9'}, ...
                {'PO5','FT9','FT10','P9','P10','PO9','PO10','IZ'}, ...
                {'PO5','AF5','AF1','T9','T10','TP9','TP10','P5'}, ...
                {'PO5','T9','T10'}, ...
                {'PO5','T9','T8','T10','TP9','TP10','I1','I2'}, ...
                {'PO5','AF7','T9','T10','TP9','TP10'}, ...
                {'PO5','T10'}, ...
                {'PO5','T9','T10','I2'}, ...
                {'PO5','T10','I2'}, ...
                {'PO5','FP1','FT9','FT10','T9','T10','PZ','I1','I2'}, ... %S20
                {'PO5','FP1','T9','T7','T10','CPZ','PO9','PO3','PO4','I2'}, ...
                {'PO5','FP1','T9','T10','CPZ','I1','I2'}, ...
                {'T10','FT9'}, ... % to a lesser extent: T7, AF7, some F's
                {'PO3','P8'},...
                {},... %25
                {'AF6','AF8'},... % usually ok, but some random bumps
                {}};
            
%             folders = {'2012-03-29-Pilot-S9', '2012-04-09-Pilot-S10', '2012-04-18-Pilot-S11', ...
%                 '2012-05-04-Pilot-S13', '2012-05-29-Pilot-S14', '2012-06-12-Pilot-S15', ...
%                 '2012-06-15-Pilot-S17', '2012-06-21-Pilot-S18', '2012-06-22-Pilot-S19',...
%                 '2012-06-25-Pilot-S20', '2013-04-15-Pilot-S21', '2013-04-16-Pilot-S22', ...
%                 '2014-02-13-Pilot-S23', '2014-02-17-Pilot-S24', '2014-02-20-Pilot-S25', ...
%                 '2014-03-04-Pilot-S26', '2014-03-10-Pilot-S27'};
% 
%             basedir = '/Users/dave/Documents/Data/Squares';

    end

    for iSubj = 1:numel(subjects)
        cd(basedir);
        cd(folders{iSubj});
        subject = subjects(iSubj);
        badElecs = badElecs_cell{iSubj};
        
        % Redo check
        if 0%exist(sprintf('%s-%d%s.set',prefix,subject,combosuffix_interp),'file')~=0
            fprintf('Skipping %s-%d interpolation...\n',prefix,subject);
            continue;
        end
        
        startfile = sprintf('%s-%d%s.set',prefix,subject,combosuffix);
        endfile = sprintf('%s-%d%s.set',prefix,subject,combosuffix_interp);
        % Remove outer electrodes
        RemoveElectrodes(startfile,endfile,outerElecs);
        % Interpolate bad electrodes & overwrite
        InterpolateElectrodes(endfile,endfile,setdiff(badElecs,outerElecs));

    end
    
end

%% Remove EOG

combosuffix_noeog = [combosuffix_interp '-noeog'];

experiment = 'Squares';
prefix = 'sq';            
subjects = [9:11 13:15 17:19 20:27];            

folders = {'2012-03-29-Pilot-S9', '2012-04-09-Pilot-S10', '2012-04-18-Pilot-S11', ...
    '2012-05-04-Pilot-S13', '2012-05-29-Pilot-S14', '2012-06-12-Pilot-S15', ...
    '2012-06-15-Pilot-S17', '2012-06-21-Pilot-S18', '2012-06-22-Pilot-S19',...
    '2012-06-25-Pilot-S20', '2013-04-15-Pilot-S21', '2013-04-16-Pilot-S22', ...
    '2014-02-13-Pilot-S23', '2014-02-17-Pilot-S24', '2014-02-20-Pilot-S25', ...
    '2014-03-04-Pilot-S26', '2014-03-10-Pilot-S27'};

basedir = '/Users/dave/Documents/Data/Squares';


for iSubj = 1:numel(subjects)
    cd(basedir);
    cd(folders{iSubj});
    subject = subjects(iSubj);
    
    % Redo check
    if exist(sprintf('%s-%d%s.set',prefix,subject,combosuffix_noeog),'file')~=0
        fprintf('Skipping %s-%d eog removal...\n',prefix,subject);
        continue;
    end

    EEG = RemoveEogComponents(prefix,subject,[],combosuffix_interp,0);
    pop_saveset(EEG,'filename',sprintf('%s-%d%s.set',prefix,subject,combosuffix_noeog));
end

%% REMOVE BLINKS FROM SQFIX DATASETS
doPlot = true;


for iExp = 1:numel(experiments)
    experiment = experiments{iExp};
    [subjects,basedir,folders,offsets] = GetSquaresSubjects(experiment);
    switch experiment
        case 'SquaresFix3'
            prefix = 'sf3';
%             subjects = 1:12;
%             offsets = zeros(1,numel(subjects));
%             basedir = '/Users/dave/Documents/Data/SquaresFix3';
%             folders = {'2013-05-15-Pilot-S1'    '2013-12-02-Pilot-S2'    '2013-12-04-Pilot-S3'    '2013-12-04-Pilot-S4',...
%                 '2013-12-05-Pilot-S5'    '2013-12-16-Pilot-S6'    '2013-12-19-Pilot-S7'    '2014-02-04-Pilot-S8',...
%                 '2014-02-06-Pilot-S9'    '2014-02-06-Pilot-S10'   '2014-02-12-Pilot-S11'   '2014-02-18-Pilot-S12'};

        case 'SquaresFix'
            prefix = 'sf';
%             subjects = [1:10 12:13];
%             offsets = zeros(1,numel(subjects));
%             basedir = '/Users/dave/Documents/Data/SquaresFix';
%             folders = {'2013-03-05-Pilot-S1','2013-03-11-Pilot-S2','2013-03-18-Pilot-S3',...
%                 '2013-03-19-Pilot-S4','2013-04-01-Pilot-S5','2013-04-02-Pilot-S6',...
%                 '2013-04-09-Pilot-S7','2013-04-11-Pilot-S8','2013-04-11-Pilot-S9',...
%                 '2013-04-17-Pilot-S10','2014-02-17-Pilot-S12','2014-03-06-Pilot-S13'};
        case 'Squares'
            continue;
    end

    for iSubj = 1:numel(subjects)
        cd(basedir);
        cd(folders{iSubj});
        subject = subjects(iSubj);
        offset_ms = offsets(iSubj);
        
        % Redo check
        if 0%exist(sprintf('%s-%d%s.set',prefix,subject,combosuffix_noeog),'file')~=0
            fprintf('Skipping %s-%d eog removal...\n',prefix,subject);
            continue;
        end
        
        % Load data        
        eegFilename = sprintf('%s-%d%s.set',prefix,subject,combosuffix_interp);
        fprintf('Loading %s...\n',eegFilename);
        EEG = pop_loadset(eegFilename);
        y = loadBehaviorData(subject,[],prefix);
        
        disp('Finding blink component...')        
        comp_blink = GetBlinkComponent(EEG,offset_ms);

        disp('Removing blink component...')
        EEG.data = SubtractOutComponents(EEG.data,comp_blink);

        if doPlot
            disp('Plotting component...')
            figure;            
            topoplot(double(comp_blink),EEG.chanlocs);
            title('Blink component')
            MakeFigureTitle(EEG.setname);
        end
        
        % Save dataset
        pop_saveset(EEG,'filename',sprintf('%s-%d%s.set',prefix,subject,combosuffix_noeog));
        disp('Done!')
    end
end

   
