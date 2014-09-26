% PlotReconstructedTrials_script.m
%
% Created 8/7/14 by DJ.

% sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','NULL','sf-Circle'};
% sequence{2} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pT_{0/2}','pT^*_{1/2}','sf-Circle'};
% sequence{3} = {'pT_{0/2}','pD_{1/2}','pD_{1/2}','pD_{1/2}','pT^*_{1/2}','sf-Circle'};
% sequence{4} = {'pT_{0/2}','pT^*_{1/2}','pD_{+/2}','pD_{+/2}','pD_{+/2}','sf-Circle'};
% sequence{5} = {'pT_{0/2}','pT^*_{1/2}','pT_{+/2}','pT_{+/2}','pT_{+/2}','sf-Circle'};
% seqnames = {'DDDDD','DDDTT','TDDDT','TTDDD','TTTTT'};

sequence{1} = {'pT_{0/2}','pD_{1/2}','pD_{1/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};
sequence{2} = {'pD_{0/2}','pT_{0/2}','pD_{1/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};
sequence{3} = {'pD_{0/2}','pD_{0/2}','pT_{0/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};
sequence{4} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pT_{0/2}','pD^*_{-/2}','sf-Circle'};
sequence{5} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','pT_{-/2}','sf-Circle'};

% sequence{1} = {'aT_{0/2}','aD_{1/2}','aD_{1/2}','aD_{1/2}','aD^*_{-/2}','Circle'};
% sequence{2} = {'aD_{0/2}','aT_{0/2}','aD_{1/2}','aD_{1/2}','aD^*_{-/2}','Circle'};
% sequence{3} = {'aD_{0/2}','aD_{0/2}','aT_{0/2}','aD_{1/2}','aD^*_{-/2}','Circle'};
% sequence{4} = {'aD_{0/2}','aD_{0/2}','aD_{0/2}','aT_{0/2}','aD^*_{-/2}','Circle'};
% sequence{5} = {'aD_{0/2}','aD_{0/2}','aD_{0/2}','aD^*_{-/2}','aT_{+/2}','Circle'};

seqnames = {'TDDDD','DTDDD','DDTDD','DDDTD','DDDDT'};

% sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
% sequence{2} = {'pD_{0/3}','pD_{0/3}','pT_{0/3}','pT_{1/3}','pT^*_{2/3}','sf3-Circle'};
% sequence{3} = {'pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pT^*_{2/3}','sf3-Circle'};
% sequence{4} = {'pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pD_{+/3}','pD_{+/3}','sf3-Circle'};
% seqnames = {'DDDDD','DDTTT','TDTDT','TTTDD'};

chansToPlot = {'FZ';'CZ';'PZ';'OZ'};
RF_seq = cell(1,length(sequence));
for i=1:length(sequence)
    [RF_seq{i},t_seq,titlestring] = PlotReconstructedTrial(group_RF_sf,tResponse,event_list_sf,sequence{i},chanlocs,chansToPlot);
end

PlotResponseFnsGrid(cat(3,RF_seq{:}),seqnames,t_seq,chanlocs,chansToPlot);
%%
for i=1:numel(chansToPlot)
    subplot(size(chansToPlot,1),size(chansToPlot,2),i);
    title('');
    ylabel(chansToPlot{i});
    ylim([-4 4]);
end

subplot(4,1,1);
title('SF3-v3.6-RampUp-ridge100, group RF')
set(gcf,'Position',[1188 354 1369 1152]);
