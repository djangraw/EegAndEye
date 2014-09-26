% Run TEMP_MakeEfficiencyPlots first!
%% Figure 1
params.doTerrain = 1;
params.doDots = 1;
params.doPath = 0;
params.doEegPTs = 0;
params.doTagPTs = 0;
params.doTagRanks = 0;
params.doTravSales = 0;
clf;
ImageAllSessions_BioNav(subject,sessions,'GridHuge.png',[15 9.5],iObjects_eeg_pt, params,usegridconstraints);
a = get(gcf,'Children');
axis(a(end),[-14.5 1207.5 -9 168])
%% Figure 2
params.doTerrain = 0;
params.doDots = 1;
params.doPath = 0;
params.doEegPTs = 0;
params.doTagPTs = 0;
params.doTagRanks = 0;
params.doTravSales = 0;
clf;
ImageAllSessions_BioNav(subject,sessions,'GridHuge.png',[15 9.5],iObjects_eeg_pt, params,usegridconstraints);
a = get(gcf,'Children');
axis(a(end),[-14.5 1207.5 -9 168])
%% Figure 3
params.doTerrain = 0;
params.doDots = 1;
params.doPath = 1;
params.doEegPTs = 0;
params.doTagPTs = 0;
params.doTagRanks = 0;
params.doTravSales = 0;
clf;
ImageAllSessions_BioNav(subject,sessions,'GridHuge.png',[15 9.5],iObjects_eeg_pt, params,usegridconstraints);
a = get(gcf,'Children');
axis(a(end),[-14.5 1207.5 -9 168])
%% Figure 4
params.doTerrain = 0;
params.doDots = 1;
params.doPath = 1;
params.doEegPTs = 1;
params.doTagPTs = 0;
params.doTagRanks = 0;
params.doTravSales = 0;
clf;
ImageAllSessions_BioNav(subject,sessions,'GridHuge.png',[15 9.5],iObjects_eeg_pt, params,usegridconstraints);
a = get(gcf,'Children');
axis(a(end),[-14.5 1207.5 -9 168])
%% Figure 5
params.doTerrain = 0;
params.doDots = 0;
params.doPath = 1;
params.doEegPTs = 0;
params.doTagPTs = 1;
params.doTagRanks = 1;
params.doTravSales = 0;
clf;
ImageAllSessions_BioNav(subject,sessions,'GridHuge.png',[15 9.5],iObjects_eeg_pt, params,usegridconstraints);
a = get(gcf,'Children');
axis(a(end),[-14.5 1207.5 -9 168])
%% Figure 6
params.doTerrain = 0;
params.doDots = 0;
params.doPath = 0;
params.doEegPTs = 0;
params.doTagPTs = 1;
params.doTagRanks = 1;
params.doTravSales = 0;
clf;
ImageAllSessions_BioNav(subject,sessions,'GridHuge.png',[15 9.5],iObjects_eeg_pt, params,usegridconstraints);
a = get(gcf,'Children');
axis(a(end),[-14.5 1207.5 -9 168])
%% Figure 7
params.doTerrain = 0;
params.doDots = 0;
params.doPath = 0;
params.doEegPTs = 0;
params.doTagPTs = 1;
params.doTagRanks = 1;
params.doTravSales = 1;
clf;
ImageAllSessions_BioNav(subject,sessions,'GridHuge.png',[15 9.5],iObjects_eeg_pt, params,usegridconstraints);
a = get(gcf,'Children');
axis(a(end),[-14.5 1207.5 -9 168])
