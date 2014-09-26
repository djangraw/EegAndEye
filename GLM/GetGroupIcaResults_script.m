% GetGroupIcaResults_script.m
%
% Created 8/13/14 by DJ.


experiment = 'sf';
subjects = [1:10 12:13];
% subjects = [9:11, 13:15 17:27];
if strcmp(experiment,'sf3')
    rules = {'T0vD0','D1vD0','T1vD0','D2vD0','T2vD0'};
    Cmap = GetSquaresEventColormap({'pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pT^*_{2/3}'});
else
    rules = {'T0vD0','D1vD0','T1vD0'};
    Cmap = GetSquaresEventColormap({'pT_{0/2}','pD_{1/2}','pT^*_{1/2}'});
end
[group_RF,group_Z,contrastFns] = deal([]);
legendstr = {};
for i=1:numel(rules)
    foo = load(sprintf('%s-v3pt5-RampUpResults-%scontrast',experiment,rules{i}));
    group_RF = cat(3,group_RF,foo.group_RF);
    group_Z = cat(3,group_Z,foo.group_Z);
    contrastFns = cat(3,contrastFns,foo.contrastFns);
    legendstr = cat(2,legendstr,foo.legendstr);
end


%% Perform ICA on group results
figure(1); clf; set(gcf,'Position',[31   627   795   879]);
[weights_grp,courses_grp,icawinv_grp,varca_grp] = ApplyIcaToGlmResults(group_RF,foo.tResponse,[0 750],foo.chanlocs,rules,6,Cmap);
MakeFigureTitle(sprintf('%s-3.5-RampuUp, ICA on group results',experiment));
%% Get & Plot ICA weights for each subject

[weights,courses,icawinv,varca] = deal(cell(1,size(contrastFns,4)));
for i=1:size(contrastFns,4)
    figure(100+i); clf;
    [weights{i},courses{i},icawinv{i},varca{i}] = ApplyIcaToGlmResults(contrastFns(:,:,:,i),foo.tResponse,[0 750],foo.chanlocs,rules,6,Cmap);
end
CascadeFigures(100+(1:size(contrastFns,4)),5);


%% Get & plot group averages
[avgW,avgTC,steW,steTC,multipliers] = MatchComponents(icawinv,courses);
figure(2); clf; set(gcf,'Position',[767   627   795   879]);
PlotGroupSvdResults(avgW,avgTC,foo.chanlocs,foo.tResponse,steW,steTC,rules,Cmap);

%% Save
save(sprintf('%s-v3pt5-RampUp-ICA',experiment),'weights*','courses*','icawinv*','varca*','avg*','ste*')

%% Get stats

% - tc is an nComponents x nTrialTypes cell matrix, each of which
% contains an nSubjects x nTimePoints matrix of doubles.

% apply multipliers
temp = cat(4,courses{:});
for i=1:size(temp,4) % subjects
    for j=1:size(temp,3) % components
        temp(:,:,j,i) = multipliers(i,j)*temp(:,:,j,i);
    end
end
% rearrange courses matrix
tc = cell(size(courses{1},3),size(courses{1},1)+1);
for i=1:size(courses{1},3) % components
    for j=1:size(courses{1},1) % trialtypes
        tc{i,j} = squeeze(temp(j,:,i,:))';
    end    
    tc{i,end} = zeros(size(tc{i,1}));
end
           
figure(3); clf; set(gcf,'Position',[1501 627 641 879]);
GetCrossSubjectStats(tc,foo.tResponse,[rules,{'zero'}],[Cmap; 0 0 0]);
