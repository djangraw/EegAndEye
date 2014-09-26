% TEMP_UpdateGroupStdErr.m
%
% Add StdErr value to saved group contrast results files.
%
% Created ~9/23/14 by DJ.

%% UPDATE GROUP STANDARD ERROR
% experiments = {'sq','sf','sf3'};
% rules = {'T0vD0','T1vT0','D1vD0','T1vD0','T+vD0','D+vD0','T*vD*','D*vD0'};
experiments = {'sf3'};
rules = {'D2vD1','D2vD0','T2vT1','T2vD0'};

suffix_in = 'Type-v3pt6-Matrices';
suffix_out = 'Type-v3pt6-10foldResults';
for iExp = 1:numel(experiments)
    experiment = experiments{iExp};
    [subjects,basedir,folders] = GetSquaresSubjects(experiment);
    cd(basedir);
    for iRule = 1:numel(rules)
        fprintf('exp %s, rule %s',experiment,rule);
        rule = rules{iRule};

        %% Get group-level stats
        load(sprintf('%s/%s-%s-%scontrast.mat',cd,experiment,suffix_out,rule));
        clear suffix
        
        [group_RF,group_P,group_SE] = RunTopLevelGlm_EEG(contrastFns,contrastVar,multcorrect);
        group_Z = norminv(group_P);
        %%
        cd(basedir)
        save(sprintf('%s-%s-%scontrast.mat',experiment,suffix_out,rule),'event_*','contrast*',...
            'group_*','chanlocs','tResponse','experiment','subjects','suffix_in','suffix_out',...
            'rule','sequence','tContrast','legendstr','lambda','multcorrect')
    end
end

%% GET COMPONENT MEAN AND STANDARD ERROR

