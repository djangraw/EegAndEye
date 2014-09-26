% RunGlm_MultiSubject_script.m
% originally called TEMP_RunGlm_MultiSubject_Full.m
%
% Created 4/13 by DJ.
% Updated 4/24/13 by DJ - band loop, comments
% Updated 3/3/14 by DJ - switched to v3.0: new sf3 subj, 40Hz files

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
        %%% ------- SQFIX3 -------- %%%
    prefix = 'sf3'; 
% %     % Type
% %     regressor_events = {{'sf3-SqNum1' 'sf3-SqNum2' 'sf3-SqNum3' 'sf3-SqNum4' 'sf-SqNum5' 'sf3-Circle'}; ...
% %         {'sf3-Dist-0T' 'sf3-Dist-1T' 'sf3-Dist-2T' 'sf3-Targ1' 'sf3-Targ2' 'sf3-Targ3' 'sf3-Incompl' 'sf3-Irrel'}}; 
% % %     suffixes = {sprintf('%sSqNum-v2pt4.mat',band_prefix); sprintf('%sSqNum-Type-v2pt4.mat',band_prefix)};
% %     suffixes = {sprintf('%sSqNum-v3pt0.mat',band_prefix); sprintf('%sSqNum-Type-v3pt0.mat',band_prefix)};
%     regressor_events = {{'sf3-Square' 'sf3-Circle'};...
%         {'sf3-SqNum1' 'sf3-SqNum2' 'sf3-SqNum3' 'sf3-SqNum4' 'sf3-SqNum5'}; ...
%         {'pD_{0/3}','pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pT^*_{2/3}','pD^*_{-/3}','pT_{+/3}','pD_{+/3}','pD_{-/3}'}}; 
%     suffixes = {sprintf('%sEvents-v3pt1.mat',band_prefix); sprintf('%sEvents-SqNum-v3pt1.mat',band_prefix); sprintf('%sEvents-SqNum-Type-v3pt1.mat',band_prefix)};
%     influence = {[0 500] [0 500] [0 750]};
%     add_regressor_events = {}; 
%     add_suffixes = {};
% 

%     regressor_events = {{'sf3-Square' 'sf3-Circle','TrialStart'};...        
%         {'pD_{0/3}','pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pT^*_{2/3}','pD^*_{-/3}','pT_{+/3}','pD_{+/3}','pD_{-/3}','pD_{0/3}-RampUp','pT_{0/3}-RampUp','pD_{1/3}-RampUp','pT_{1/3}-RampUp','pD_{2/3}-RampUp','pT^*_{2/3}-RampUp'}}; 
%     suffixes = {sprintf('%sEvents-v3pt4-RampUp.mat',band_prefix); sprintf('%sEvents-Type-v3pt4-RampUp.mat',band_prefix)};    
%     influence = {[0 750] [0 750]};
%     add_regressor_events = {}; 
%     add_suffixes = {};

%     regressor_events = {{'sf3-Square' 'sf3-Circle','TrialStart','pD_{0/3}','pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pT^*_{2/3}','pD^*_{-/3}','pT_{+/3}','pD_{+/3}','pD_{-/3}','pT_{-/3}'...
%         'pD_{0/3}-RampUp','pT_{0/3}-RampUp','pD_{1/3}-RampUp','pT_{1/3}-RampUp','pD_{2/3}-RampUp','pT^*_{2/3}-RampUp','pD^*_{-/3}-RampUp','pT_{+/3}-RampUp','pD_{+/3}-RampUp','pD_{-/3}-RampUp','pT_{-/3}-RampUp'}}; 
%     suffixes = {sprintf('%sType-v3pt6-RampUp.mat',band_prefix)};    
%     influence = {[0 1000]};
%     add_regressor_events = {}; 
%     add_suffixes = {};
    
    regressor_events = {{'sf3-Square' 'sf3-Circle','TrialStart','pD_{0/3}','pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pT^*_{2/3}','pD^*_{-/3}','pT_{+/3}','pD_{+/3}','pD_{-/3}','pT_{-/3}'...
    'pD_{0/3}-Peak','pT_{0/3}-Peak','pD_{1/3}-Peak','pT_{1/3}-Peak','pD_{2/3}-Peak','pT^*_{2/3}-Peak','pD^*_{-/3}-Peak','pT_{+/3}-Peak','pD_{+/3}-Peak','pD_{-/3}-Peak','pT_{-/3}-Peak'}}; 
    suffixes = {sprintf('%sType-v3pt6-Peak.mat',band_prefix)};    
    influence = {[0 1000]};
    add_regressor_events = {}; 
    add_suffixes = {};



%     % SqNum
% %     regressor_events = {{'sf3-Square' 'sf3-Circle'};...
% %         {'pD_{0/3}','pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pT^*_{2/3}','pD^*_{-/3}','pT_{+/3}','pD_{+/3}','pD_{-/3}'}; ...
% %         {'sf-SqNum1' 'sf-SqNum2' 'sf-SqNum3' 'sf-SqNum4' 'sf-SqNum5'}}; 
% %     suffixes = {sprintf('%sEvents-v3pt1.mat',band_prefix); sprintf('%sEvents-Type-v3pt1.mat',band_prefix); sprintf('%sEvents-Type-SqNum-v3pt1.mat',band_prefix)};
% %     influence = {[0 500] [0 500] [0 500]};
% %     add_regressor_events = {}; 
% %     add_suffixes = {};
    %%% ------- ----- -------- %%%

    
    
    
    %%% ------- SQFIX -------- %%%
%     prefix = 'sf'; 
%     % Type
% %     regressor_events = {{'sf-SqNum1' 'sf-SqNum2' 'sf-SqNum3' 'sf-SqNum4' 'sf-SqNum5' 'sf-Circle'}; ...
% %         {'sf-Dist-0T' 'sf-Dist-1T' 'sf-Integ' 'sf-Compl' 'sf-Incompl' 'sf-Irrel'}}; 
% %     suffixes = {sprintf('%sSqNum-v2pt4.mat',band_prefix); sprintf('%sSqNum-Type-v2pt4.mat',band_prefix)};
% %     suffixes = {sprintf('%sSqNum-v3pt0.mat',band_prefix); sprintf('%sSqNum-Type-v3pt0.mat',band_prefix)};
% 
% %     regressor_events = {{'sf-Square' 'sf-Circle'};...
% %         {'sf-SqNum1' 'sf-SqNum2' 'sf-SqNum3' 'sf-SqNum4' 'sf-SqNum5'}; ...
% %         {'pD_{0/2}' 'pT_{0/2}' 'pD_{1/2}', 'pT^*_{1/2}' 'pD^*_{-/2}' 'pT_{+/2}' 'pD_{+/2}'}}; 
% %     suffixes = {sprintf('%sEvents-v3pt1.mat',band_prefix); sprintf('%sEvents-SqNum-v3pt1.mat',band_prefix); sprintf('%sEvents-SqNum-Type-v3pt1.mat',band_prefix)};
%     
% %     regressor_events = {{'sf-Square' 'sf-Circle','TrialStart'};...        
% %         {'pD_{0/2}' 'pT_{0/2}' 'pD_{1/2}', 'pT^*_{1/2}' 'pD^*_{-/2}' 'pT_{+/2}' 'pD_{+/2}','pD_{0/2}-RampUp','pT_{0/2}-RampUp','pD_{1/2}-RampUp','pT^*_{1/2}-RampUp'}}; 
% %     suffixes = {sprintf('%sEvents-v3pt4-RampUp.mat',band_prefix); sprintf('%sEvents-Type-v3pt4-RampUp.mat',band_prefix)};
% %     influence = {[0 750] [0 750]};
% %     regressor_events = {{'sf-Square','sf-Circle','TrialStart','pD_{0/2}' 'pT_{0/2}' 'pD_{1/2}', 'pT^*_{1/2}' 'pD^*_{-/2}' 'pT_{+/2}' 'pD_{+/2}','pD_{0/2}-RampUp','pT_{0/2}-RampUp','pD_{1/2}-RampUp','pT^*_{1/2}-RampUp'}}; 
% %     suffixes = {sprintf('%sType-v3pt5-RampUp.mat',band_prefix)};
% %     influence = {[0 750]};
%     regressor_events = {{'sf-Square','sf-Circle','TrialStart','pD_{0/2}' 'pT_{0/2}' 'pD_{1/2}', 'pT^*_{1/2}' 'pD^*_{-/2}' 'pT_{+/2}' 'pD_{+/2}','pD_{-/2}','pT_{-/2}',...
%         'pD_{0/2}-Peak','pT_{0/2}-Peak','pD_{1/2}-Peak','pT^*_{1/2}-Peak','pD^*_{-/2}-Peak' 'pT_{+/2}-Peak' 'pD_{+/2}-Peak'}}; 
%     suffixes = {sprintf('%sType-v3pt6-Peak.mat',band_prefix)};
%     influence = {[0 1000]};
%     add_regressor_events = {}; 
%     add_suffixes = {};

%     % SqNum
%     regressor_events = {{'sf-Square' 'sf-Circle'};...
%         {'pD_{0/2}' 'pT_{0/2}' 'pD_{1/2}', 'pT^*_{1/2}' 'pD^*_{-/2}' 'pT_{+/2}' 'pD_{+/2}'}; ...
%         {'sf-SqNum1' 'sf-SqNum2' 'sf-SqNum3' 'sf-SqNum4' 'sf-SqNum5'}}; 
%     suffixes = {sprintf('%sEvents-v3pt1.mat',band_prefix); sprintf('%sEvents-Type-v3pt1.mat',band_prefix); sprintf('%sEvents-Type-SqNum-v3pt1.mat',band_prefix)};
%     influence = {[0 500] [0 500] [0 500]};
%     add_regressor_events = {}; 
%     add_suffixes = {};
    %%% ------- ----- -------- %%%

    

    %%% ------- SQUARES -------- %%%
%     prefix = 'sq';
% % %     % Type
% % % %     regressor_events = {{'LSaccade' 'RSaccade' 'TrialStart-L' 'TrialStart-R' 'TrialEnd-L' 'TrialEnd-R'};...
% % % %         {'SqNum1' 'SqNum2' 'SqNum3' 'SqNum4' 'SqNum5' 'Circle'}; ...
% % % %         {'New-Dist-0T' 'New-Dist-1T' 'New-Integ' 'Compl' 'Incompl' 'New-Irrel'}}; 
% % % %     suffixes = {sprintf('%sLREvents-v2pt4%s.mat',band_prefix,contrast_suffix); ...
% % % %         sprintf('%sLREvents-SqNum-v2pt4%s.mat',band_prefix,contrast_suffix); ...
% % % %         sprintf('%sLREvents-SqNum-Type-v2pt4%s.mat',band_prefix,contrast_suffix)};
% % %     regressor_events = {{'Square' 'Circle' 'TrialStart' 'TrialEnd'};...
% % %         {'SqNum1' 'SqNum2' 'SqNum3' 'SqNum4' 'SqNum5'}; ...
% % %         {'aD_{0/2}' 'aT_{0/2}' 'aD_{1/2}', 'aT^*_{1/2}' 'aD^*_{-/2}' 'aT_{+/2}' 'aD_{+/2}'}}; 
% % %     suffixes = {sprintf('%sEvents-v3pt1%s.mat',band_prefix,contrast_suffix); ...
% % %         sprintf('%sEvents-SqNum-v3pt1%s.mat',band_prefix,contrast_suffix); ...
% % %         sprintf('%sEvents-SqNum-Type-v3pt1%s.mat',band_prefix,contrast_suffix)};
% % %     influence = {[0 750] [0 750] [0 750]};
% % %     add_regressor_events = {}; 
% % %     add_suffixes = {};
% 
% %     regressor_events = {{'Square' 'Circle','TrialStart'};...        
% %         {'aD_{0/2}' 'aT_{0/2}' 'aD_{1/2}', 'aT^*_{1/2}' 'aD^*_{-/2}' 'aT_{+/2}' 'aD_{+/2}','aD_{0/2}-RampUp','aT_{0/2}-RampUp','aD_{1/2}-RampUp','aT^*_{1/2}-RampUp'}}; 
% %     suffixes = {sprintf('%sEvents-v3pt4-RampUp.mat',band_prefix); sprintf('%sEvents-Type-v3pt4-RampUp.mat',band_prefix)};
% %     regressor_events = {{'Circle','TrialStart','aD_{0/2}' 'aT_{0/2}' 'aD_{1/2}', 'aT^*_{1/2}' 'aD^*_{-/2}' 'aT_{+/2}' 'aD_{+/2}','aD_{0/2}-RampUp','aT_{0/2}-RampUp','aD_{1/2}-RampUp','aT^*_{1/2}-RampUp'}}; 
% %     suffixes = {sprintf('%sType-v3pt5-RampUp.mat',band_prefix)};
% %     influence = {[0 750]};
% %     add_regressor_events = {}; 
% %     add_suffixes = {};
%     regressor_events = {{'Square','Circle','TrialStart','aD_{0/2}' 'aT_{0/2}' 'aD_{1/2}', 'aT^*_{1/2}' 'aD^*_{-/2}' 'aT_{+/2}' 'aD_{+/2}','aD_{-/2}','aT_{-/2}',...
%         'aD_{0/2}-Peak','aT_{0/2}-Peak','aD_{1/2}-Peak','aT^*_{1/2}-Peak','aD^*_{-/2}-Peak' 'aT_{+/2}-Peak' 'aD_{+/2}-Peak'}}; 
%     suffixes = {sprintf('%sType-v3pt6-Peak.mat',band_prefix)};
%     influence = {[0 1000]};
%     add_regressor_events = {}; 
%     add_suffixes = {};


%   
%     % SqNum
%     regressor_events = {{'Square' 'Circle' 'TrialStart' 'TrialEnd'};...
%         {'aD_{0/2}' 'aT_{0/2}' 'aD_{1/2}', 'aT^*_{1/2}' 'aD^*_{-/2}' 'aT_{+/2}' 'aD_{+/2}'}; ...
%         {'SqNum1' 'SqNum2' 'SqNum3' 'SqNum4' 'SqNum5'}}; 
%     suffixes = {sprintf('%sEvents-v3pt1%s.mat',band_prefix,contrast_suffix); ...
%         sprintf('%sEvents-Type-v3pt1%s.mat',band_prefix,contrast_suffix); ...
%         sprintf('%sEvents-Type-SqNum-v3pt1%s.mat',band_prefix,contrast_suffix)};
%     influence = {[0 750] [0 750] [0 750]};
% %     add_regressor_events = {}; 
% %     add_suffixes = {};
    %%% ------- ----------- ------- %%%

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

    switch params.prefix
        case 'sq'
            basedir = '/Users/dave/Documents/Data/Squares';
            subjects = [9:11, 13:15, 17:19 20:27];
%             offsets = [42 102 -12 -18 -12 -12 6 -12 18 0 0 0 0 0 0 0 0];
            offsets = [30 100 -10 -30 -10 -30 0 0 0 -10 0 40 -20 -20 -10  10  20];
            folders = {'2012-03-29-Pilot-S9', '2012-04-09-Pilot-S10', '2012-04-18-Pilot-S11', ...
                '2012-05-04-Pilot-S13', '2012-05-29-Pilot-S14', '2012-06-12-Pilot-S15', ...
                '2012-06-15-Pilot-S17', '2012-06-21-Pilot-S18', '2012-06-22-Pilot-S19',...
                '2012-06-25-Pilot-S20', '2013-04-15-Pilot-S21', '2013-04-16-Pilot-S22', ...
                '2014-02-13-Pilot-S23', '2014-02-17-Pilot-S24', '2014-02-20-Pilot-S25', ...
                '2014-03-04-Pilot-S26', '2014-03-10-Pilot-S27'};
            % Use saccade-based rejection rules
            params.trial_rej_rules = {'skipped_ends'  'skip'  'backward' 'early_button' 'late_button' 'wrong_button'}; % added button stuff for v2.3
            params.artifact_events = {'Button'  'BlinkEnd'  'BlinkStart'  'Errant'  'Cross'};
            params.influence = influence;%[0 750];
            
        case 'sf'
            basedir = '/Users/dave/Documents/Data/SquaresFix';
            subjects = [1:10 12:13];
%             offsets = zeros(1,numel(subjects));
            offsets = [40 30 70 70 20 0 20 -20 -20 -10 -10 0];
            folders = {'2013-03-05-Pilot-S1','2013-03-11-Pilot-S2','2013-03-18-Pilot-S3',...
                '2013-03-19-Pilot-S4','2013-04-01-Pilot-S5','2013-04-02-Pilot-S6',...
                '2013-04-09-Pilot-S7','2013-04-11-Pilot-S8','2013-04-11-Pilot-S9',...
                '2013-04-17-Pilot-S10','2014-02-17-Pilot-S12','2014-03-06-Pilot-S13'};
            % Only use non-saccade-based rejection rules
            params.trial_rej_rules = {'early_button' 'late_button' 'wrong_button'}; % added button stuff for v2.3
            params.artifact_events = {'Button'  'BlinkEnd'  'BlinkStart'};
            params.influence = influence;%[0 500];
        case 'sf3'
            basedir = '/Users/dave/Documents/Data/SquaresFix3';
            subjects = 1:12;
%             offsets = zeros(1,numel(subjects));
            offsets = [-20 50 90 80 80 -10 0 -30 20 0 30 20]; % added 7/21/14
            folders = {'2013-05-15-Pilot-S1'    '2013-12-02-Pilot-S2'    '2013-12-04-Pilot-S3'    '2013-12-04-Pilot-S4',...
                '2013-12-05-Pilot-S5'    '2013-12-16-Pilot-S6'    '2013-12-19-Pilot-S7'    '2014-02-04-Pilot-S8',...
                '2014-02-06-Pilot-S9'    '2014-02-06-Pilot-S10'   '2014-02-12-Pilot-S11'   '2014-02-18-Pilot-S12'};
            % Only use non-saccade-based rejection rules
            params.trial_rej_rules = {'early_button' 'late_button' 'wrong_button'}; % added button stuff for v2.3
            params.artifact_events = {'Button'  'BlinkEnd'  'BlinkStart'};
            params.influence = influence;%[0 500];
    end
    
    % --- RUN IT
    cd(basedir)
    % band_suffix = '-lower1alpha';
    % band_suffix = '-theta';
    % band_suffix = '';

    for i=1:numel(subjects)
        fprintf('--- SUBJECT %d ---\n',subjects(i));
        cd(folders{i});
        subject = subjects(i);
        params.offset = offsets(i);
        % Run Type GLM
%         if strcmp(sprintf('%s-%d-GLMresults-%s',params.prefix,subject,suffixes{end}),'sq-17-GLMresults-UpperAlpha_v2-LREvents-SqNum-Type-v2pt4IncvsD1.mat')
%             disp('***SKIP THIS ONE!***')
%             cd ..
%             continue;
%         end
        if ~exist(sprintf('%s-%d-GLMresults-%s',params.prefix,subject,suffixes{end}),'file')
%             RunGlmLikeGui(sprintf('%s-%d-all-filtered-50Hz-interpduds%s.set',params.prefix,subject,band_suffix),regressor_events,suffixes,params)
%             if strcmp(params.prefix,'sq')
                RunGlmLikeGui(sprintf('%s-%d-all-40Hz-fs100-interp-noeog%s.set',params.prefix,subject,band_suffix),regressor_events,suffixes,params)
%             else
%                 RunGlmLikeGui(sprintf('%s-%d-all-40Hz-fs100-interp%s.set',params.prefix,subject,band_suffix),regressor_events,suffixes,params)
%             end
        else
            disp('--- SKIPPING REGULAR GLM ---');
        end
        % Run SqNum GLM (reusing 1st level of NewType)
        if ~isempty(add_suffixes) && ~exist(sprintf('%s-%d-GLMresults-%s',params.prefix,subject,add_suffixes{end}),'file');
%             RunAddOnLevelOfGlm(sprintf('%s-%d-GLMresults-%s',params.prefi
%             x,subject,add_suffixes{1}),1,add_regressor_events(2:end), add_suffixes(2:end));
            RunAddOnLevelOfGlm(sprintf('%s-%d-GLMresults-%s',params.prefix,subject,suffixes{end}),1,add_regressor_events(2:end), add_suffixes(2:end));
        else
            disp('--- SKIPPING ADD-ON GLM ---');        
        end
        cd ..
        fprintf('--- DONE! ---\n');
    end

end % band

end % contrast