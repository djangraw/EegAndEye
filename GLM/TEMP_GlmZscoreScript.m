% TEMP_GlmZscoreScript.m
%
% Created 6/27/12 by DJ for one-time use

S.subject = 13;
% S.rule = 'SqNum-Type-v1pt5';
S.rule = 'Type-SqNum-v1pt5';
switch S.rule
    case 'SqNum-Type-v1pt5';
        S.duringrule = 'SqNum-v1pt5';
    case 'Type-SqNum-v1pt5'
        S.duringrule = 'Type-v1pt5';
end

%%

if strcmp(S.EEG.filename,sprintf('sq-%d-GLMstart.set',S.subject))
    disp('Using loaded GMLstart set...')
elseif exist(sprintf('sq-%d-GLMstart.set',S.subject),'file')
    disp('Loading GLMstart set...')
    S.EEG = pop_loadset(sprintf('sq-%d-GLMstart.set',S.subject));
else
    disp('Making GLMstart set...')
%     EEG = pop_loadset(sprintf('sq-%d-all-filtered-50Hz-noduds.set',S.subject));
    EEG = pop_loadset(sprintf('sq-%d-all-noeog-50Hz-noduds.set',S.subject));
    EEG = SetUpGlm(S.subject,EEG,S.rule);
    EEG = pop_saveset(EEG,sprintf('sq-%d-GLMstart.set',S.subject));
    S.EEG = EEG;
end
    
disp('Loading GLM results...')
S.during = load(sprintf('sq-%d-GLMresults-%s',S.subject,S.duringrule));
S.after = load(sprintf('sq-%d-GLMresults-%s',S.subject,S.rule));
disp('Done!');

%%
PlotSquaresStats(S.EEG,S.during,S.after,S.rule,1:5);