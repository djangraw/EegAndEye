% [weights, courses] = ApplyIcaToGlmResults(cat(3,group_RF{:}),R(1).tResponse,[0 500],R(1).EEG.chanlocs,rowlabels,5);



R = R_sf3_type;
events = R(1).regressor_events{end}(1:6);

%% Run Hierarchical GLM
[contrastFns, group_RF, contrastZ, group_z] = deal(cell(1,numel(plus_events)));
for i=1:numel(events)
    [contrastFns{i}, group_RF{i}, contrastZ{i}, group_z{i}] = RunTopLevelGlm_Full(R,events{i},'',baseline_win,multcorrect);  
end
%% Apply ICA to group results
nComps = 5; % # of components to plot
[weights, courses, icawinv, varca] = ApplyIcaToGlmResults(cat(3,group_RF{:}),R(1).tResponse,[0 500],R(1).EEG.chanlocs,events,nComps);
%% Re-plot components
clf;
Cmap = GetSquaresEventColormap(events);

PlotSvdWeightsAndCourses(cat(3,group_RF{:}),R(1).tResponse,weights,varca,R(1).EEG.chanlocs,events,nComps,icawinv,Cmap);

%% Get individual components
nComps = 5; % # of components to plot
Cmap = GetSquaresEventColormap(events);
[weights, courses, icawinv, varca] = deal(cell(1,numel(R)));
for i=1:numel(R)
    iEvents = find(ismember(R(i).regressor_events{end},events));
    RF = R(i).responseFns{end}(:,:,iEvents);
    figure(200+i)
    [weights{i}, courses{i}, icawinv{i}, varca{i}] = ApplyIcaToGlmResults(RF,R(i).tResponse,[0 500],R(i).EEG.chanlocs,events,nComps,Cmap);  
    
end

%% Flip signs as necessary and plot results
[avgW, avgTC, steW, steTC] = MatchComponents(icawinv,courses);

figure;
PlotGroupSvdResults(avgW, avgTC, R(1).EEG.chanlocs, R(1).tResponse, steW, steTC, events,Cmap)

    
    