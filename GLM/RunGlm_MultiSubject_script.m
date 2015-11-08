% RunGlm_MultiSubject_script.m
% originally called TEMP_RunGlm_MultiSubject_Full.m
%
% Created 4/13 by DJ.
% Updated 4/24/13 by DJ - band loop, comments
% Updated 3/3/14 by DJ - switched to v3.0: new sf3 subj, 40Hz files
% Updated 1/27/15 by DJ - use GetSquaresSubjects(), cleaned up

% Set up Analysis
% bandnames = {'', 'Theta', 'Lower1Alpha', 'Lower2Alpha', 'UpperAlpha'}; % capitalize first letter
% bandnames = {'Beta', 'LowerGamma'}; % capitalize first letter
% bandnames = {'','Theta_v2','UpperAlpha_v2','Lower1Alpha_v2','Lower2Alpha_v2','Beta_v2','LowerGamma_v2'};
bandnames = {''};

% contrast_suffix_cell = {'CvsD1' 'CvsI','IvsD0','CvsD0','CvsInc','IncvsD1','D1vsD0'};
% contrast_events_cell = {{'Compl','New-Dist-1T'},{'Compl','New-Integ'},...
%     {'New-Integ','New-Dist-0T'},{'Compl','New-Dist-0T'},...
%     {'Compl','Incompl'},{'Incompl','New-Dist-1T'},{'New-Dist-1T','New-Dist-0T'}};

% contrast_suffix_cell = {'T1vsD1' 'T1vsT0',...
%     'T0vsD0','T1vsD0',...
%     'T1vsD-','D-vsD1',...
%     'D1vsD0','T+vsT1','T+vsT0'};
% contrast_events_cell = {{'aT^*_{1/2}','aD_{1/2}'},{'aT^*_{1/2}','aT_{0/2}'},...
%     {'aT_{0/2}','aD_{0/2}'},{'aT^*_{1/2}','aD_{0/2}'},...
%     {'aT^*_{1/2}','aD^*_{-/2}'},{'aD^*_{-/2}','aD_{1/2}'},...
%     {'aD_{1/2}','aD_{0/2}'},{'aT_{+/2}','aT^*_{1/2}'},{'aT_{+/2}','aT_{0/2}'}};

contrast_suffix_cell = {''};
contrast_events_cell = {{}};

for iBand = 1:numel(bandnames)
fprintf('--- Band: %s ---\n',bandnames{iBand});

for iContrast = 1:numel(contrast_suffix_cell)    
    contrast_suffix = contrast_suffix_cell{iContrast};
    contrast_events = contrast_events_cell{iContrast};    
    fprintf('--- Contrast: %s ---\n',contrast_suffix);
    

    
    % Set up

    band_name = bandnames{iBand}; 
    if isempty(band_name)
        [band_suffix, band_prefix] = deal('');
    else
        band_suffix = sprintf('-%s',lower(band_name));
        band_prefix = sprintf('%s-',band_name);
    end
    % DECLARE EXPERIMENT AND MODEL!!!
    prefix = 'sf';
    model = 'Type';        
     
    switch prefix        
        case 'sf3'
        %%% ------- SQFIX3 -------- %%%
            switch model
                case 'Type'
                    % 1/27/15: latest TYPE
                    regressor_events = {{'sf3-Square' 'sf3-Circle','TrialStart','pD_{0/3}','pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pT^*_{2/3}','pD^*_{-/3}','pT_{+/3}','pD_{+/3}','pD_{-/3}','pT_{-/3}'...
                        'pD_{0/3}-Peak','pT_{0/3}-Peak','pD_{1/3}-Peak','pT_{1/3}-Peak','pD_{2/3}-Peak','pT^*_{2/3}-Peak','pD^*_{-/3}-Peak','pT_{+/3}-Peak','pD_{+/3}-Peak','pD_{-/3}-Peak','pT_{-/3}-Peak'}}; 
                    suffixes = {sprintf('%sType-v3pt6-Peak.mat',band_prefix)};
                case 'SqNum'
                    % 1/27/15: new SQNUM
                    regressor_events = {{'sf3-Circle','TrialStart','sf3-SqNum1' 'sf3-SqNum2' 'sf3-SqNum3' 'sf3-SqNum4' 'sf3-SqNum5'}}; 
                    suffixes = {sprintf('%sSqNum-v3pt6.mat',band_prefix)};                    
                case 'TargDis'
                    % 1/28/15: new TARGDIS
                    regressor_events = {{'sf-Circle','TrialStart','pD_{x/3}' 'pT_{x/3}'}}; 
                    suffixes = {sprintf('%sTargDis-v3pt6.mat',band_prefix)};                    
            end
            influence = {[0 1000]};

    add_regressor_events = {}; 
    add_suffixes = {};
    %%% ------- ----- -------- %%%

    
    
        case 'sf'
        %%% ------- SQFIX -------- %%%
            switch model
                case 'Type'
                    % 1/27/15: latest TYPE
                    regressor_events = {{'sf-Square','sf-Circle','TrialStart','pD_{0/2}' 'pT_{0/2}' 'pD_{1/2}', 'pT^*_{1/2}' 'pD^*_{-/2}' 'pT_{+/2}' 'pD_{+/2}','pD_{-/2}','pT_{-/2}'}}; 
                    suffixes = {sprintf('%sType-v3pt6.mat',band_prefix)};
                case 'Peak'
                    % 1/27/15: latest TYPE
                    regressor_events = {{'sf-Square','sf-Circle','TrialStart','pD_{0/2}' 'pT_{0/2}' 'pD_{1/2}', 'pT^*_{1/2}' 'pD^*_{-/2}' 'pT_{+/2}' 'pD_{+/2}','pD_{-/2}','pT_{-/2}',...
                        'pD_{0/2}-Peak','pT_{0/2}-Peak','pD_{1/2}-Peak','pT^*_{1/2}-Peak','pD^*_{-/2}-Peak' 'pT_{+/2}-Peak' 'pD_{+/2}-Peak'}}; 
                    suffixes = {sprintf('%sType-v3pt6-Peak.mat',band_prefix)};
                case 'SqNum'
                    % 1/27/15: new SQNUM
                    regressor_events = {{'sf-Circle','TrialStart','sf-SqNum1' 'sf-SqNum2' 'sf-SqNum3' 'sf-SqNum4' 'sf-SqNum5'}}; 
                    suffixes = {sprintf('%sSqNum-v3pt6.mat',band_prefix)};
                case 'TargDis'
                    % 1/28/15: new TARGDIS
                    regressor_events = {{'sf-Circle','TrialStart','pD_{x/2}' 'pT_{x/2}'}}; 
                    suffixes = {sprintf('%sTargDis-v3pt6.mat',band_prefix)};
            end
            influence = {[0 1000]};

    add_regressor_events = {}; 
    add_suffixes = {};
    %%% ------- ----- -------- %%%

        case 'sq'
        %%% ------- SQUARES -------- %%%

            switch model
                case 'Peak'
                    % 1/27/15: latest TYPE-PEAK
                    regressor_events = {{'Square','Circle','TrialStart','aD_{0/2}' 'aT_{0/2}' 'aD_{1/2}', 'aT^*_{1/2}' 'aD^*_{-/2}' 'aT_{+/2}' 'aD_{+/2}','aD_{-/2}','aT_{-/2}',...
                        'aD_{0/2}-Peak','aT_{0/2}-Peak','aD_{1/2}-Peak','aT^*_{1/2}-Peak','aD^*_{-/2}-Peak' 'aT_{+/2}-Peak' 'aD_{+/2}-Peak'}}; 
                    suffixes = {sprintf('%sType-v3pt6-Peak.mat',band_prefix)};
                case 'SqNum'
                    % 1/27/15: new SQNUM
                    regressor_events = {{'Circle','TrialStart','SqNum1' 'SqNum2' 'SqNum3' 'SqNum4' 'SqNum5'}}; 
                    suffixes = {sprintf('%sSqNum-v3pt6.mat',band_prefix)};         
                case 'TargDis'
                    % 1/28/15: new TARGDIS
                    regressor_events = {{'sf-Circle','TrialStart','aD_{x/2}' 'aT_{x/2}'}}; 
                    suffixes = {sprintf('%sTargDis-v3pt6.mat',band_prefix)};
            end
            influence = {[0 1000]};

    add_regressor_events = {}; 
    add_suffixes = {};
    %%% ------- ----------- ------- %%%
    end
    
    
    % params = struct('offset',0,'influence',500,'stddev',0,'vthresh',75,'method','mvregress','trial_rej_rules',{'skipped_ends'  'skip'  'backward'},'artifact_events',{'Button'  'BlinkEnd'  'BlinkStart'  'Errant'  'Cross'});
    params.offset=0;
%     params.influence=[0 500]; % [0 500] for v2.3
    params.artifact_influence = [-500 500]; % added for v2.3
    params.stddev=0;
    params.vthresh=75;
    params.method='leastsquares';%'mvregress';
%     params.method='ridge';
%     params.lambda = {'auto'};
    params.trial_rej_rules = {'skipped_ends'  'skip'  'backward' 'early_button' 'late_button' 'wrong_button'}; % added button stuff for v2.3
    params.artifact_events = {'Button'  'BlinkEnd'  'BlinkStart'  'Errant'  'Cross'};
    params.prefix = prefix;    
    params.contrast_events = contrast_events;

    [subjects,basedir,folders,offsets] = GetSquaresSubjects(prefix);
    switch params.prefix
        case 'sq'
            % Use saccade-based rejection rules
            params.trial_rej_rules = {'skipped_ends'  'skip'  'backward' 'early_button' 'late_button' 'wrong_button'}; % added button stuff for v2.3
            params.artifact_events = {'Button'  'BlinkEnd'  'BlinkStart'  'Errant'  'Cross'};
            params.influence = influence;%[0 750];
        case 'sf'
            % Only use non-saccade-based rejection rules
            params.trial_rej_rules = {'early_button' 'late_button' 'wrong_button'}; % added button stuff for v2.3
            params.artifact_events = {'Button'  'BlinkEnd'  'BlinkStart'};
            params.influence = influence;%[0 500];
        case 'sf3'
            % Only use non-saccade-based rejection rules
            params.trial_rej_rules = {'early_button' 'late_button' 'wrong_button'}; % added button stuff for v2.3
            params.artifact_events = {'Button'  'BlinkEnd'  'BlinkStart'};
            params.influence = influence;%[0 500];
    end
    
    % --- RUN IT
    cd(basedir)    

    for i=1:numel(subjects)
        fprintf('--- SUBJECT %d ---\n',subjects(i));
        cd(folders{i});
        subject = subjects(i);
        params.offset = offsets(i);
        % Run Type GLM
        if ~exist(sprintf('%s-%d-GLMresults-%s',params.prefix,subject,suffixes{end}),'file')
            RunGlmLikeGui(sprintf('%s-%d-all-40Hz-fs100-interp-noeog%s.set',params.prefix,subject,band_suffix),regressor_events,suffixes,params)
        else
            disp('--- SKIPPING REGULAR GLM ---');
        end
        % Run 2nd-level GLM (reusing 1st level already saved out)
        if ~isempty(add_suffixes) && ~exist(sprintf('%s-%d-GLMresults-%s',params.prefix,subject,add_suffixes{end}),'file');
            RunAddOnLevelOfGlm(sprintf('%s-%d-GLMresults-%s',params.prefix,subject,suffixes{end}),1,add_regressor_events(2:end), add_suffixes(2:end));
        else
            disp('--- SKIPPING ADD-ON GLM ---');        
        end
        cd ..
        fprintf('--- DONE! ---\n');
    end

end % band

end % contrast