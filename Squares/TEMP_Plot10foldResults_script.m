% TEMP_Plot10foldResults_script.m
%
% Created ~9/18/14 by DJ.

[group_RF,group_Z,group_P,contrastFns] = deal([]);
legendstr = {};
experiments = {'sq','sf','sf3'};
rule = 'T1vT0';
for i=1:numel(experiments)    
    cd(basedirs{i});
    foo = load(sprintf('%s-Type-v3pt6-RampUp-10foldResults-%scontrast',experiments{i},rule));
    group_RF = cat(3,group_RF,foo.group_RF);
    group_Z = cat(3,group_Z,foo.group_Z);
    group_P = cat(3,group_P,foo.group_P);
%     contrastFns = cat(3,contrastFns,foo.contrastFns);
%     legendstr = cat(2,legendstr,foo.legendstr);
    legendstr = cat(2,legendstr,sprintf('%s, %s',experiments{i},rule));
end

%%
[group_RF,group_Z,group_P,contrastFns] = deal([]);
legendstr = {};
experiment = 'sf3';
iExp = find(strcmp(experiment,experiments));

% rules = {'D1vD0','T0vD0','T1vD0','T+vD0'};
rules = {'D2vD0','D1vD0','T0vD0','T1vD0','T2vD0','T+vD0'};
for i=1:numel(rules)    
    cd(basedirs{iExp});
    foo = load(sprintf('%s-Type-v3pt6-RampUp-10foldResults-%scontrast',experiment,rules{i}));
    group_RF = cat(3,group_RF,foo.group_RF);
    group_Z = cat(3,group_Z,foo.group_Z);
    group_P = cat(3,group_P,foo.group_P);
%     contrastFns = cat(3,contrastFns,foo.contrastFns);
%     legendstr = cat(2,legendstr,foo.legendstr);
    legendstr = cat(2,legendstr,sprintf('%s, %s',experiment,rules{i}));
end