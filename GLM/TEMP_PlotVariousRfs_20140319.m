%% Load and set up
% load('TEMP_GLMresults_2014-03-20.mat','*_type');
chansToPlot = {'F3','FZ','F4';'C3','CZ','C4';'P3','PZ','P4'};
iLevel_vec = [3 3 3];
addSquareRf = false;
%% First 4
PlotResponseFnsForThesis({R_sq_type,R_sf_type,R_sf3_type},{1:4,1:4,1:4},chansToPlot,iLevel_vec,addSquareRf);
%% Up Until Decision
PlotResponseFnsForThesis({R_sq_type,R_sf_type,R_sf3_type},{1:4,1:4,1:6},chansToPlot,iLevel_vec,addSquareRf);
%% Around Target
PlotResponseFnsForThesis({R_sq_type,R_sf_type,R_sf3_type},{3:5,3:5,5:7},chansToPlot,iLevel_vec,addSquareRf);
%% Distractors over time
PlotResponseFnsForThesis({R_sq_type,R_sf_type,R_sf3_type},{1:2:7,1:2:7,1:2:9},chansToPlot,iLevel_vec,addSquareRf);
%% Targets over time
PlotResponseFnsForThesis({R_sq_type,R_sf_type,R_sf3_type},{2:2:6,2:2:6,2:2:8},chansToPlot,iLevel_vec,addSquareRf);
%% Extra Targets
PlotResponseFnsForThesis({R_sq_type,R_sf_type,R_sf3_type},{[1:4 6], [1:4 6], [1:6 8]},chansToPlot,iLevel_vec,addSquareRf);
%% Decision Points
PlotResponseFnsForThesis({R_sq_type,R_sf_type,R_sf3_type},{2:7, 2:7, 4:9},chansToPlot,iLevel_vec,addSquareRf);
%% SqNum
% PlotResponseFnsForThesis({R_sq_type,R_sf_type,R_sf3_type},{1:5,1:5,1:5},[],[2,2,2]);
PlotResponseFnsForThesis({R_sq_sqnum,R_sf_sqnum,R_sf3_sqnum},{1:5,1:5,1:5},chansToPlot,iLevel_vec,addSquareRf);



