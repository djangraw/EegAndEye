function [contrasts, legendstr, sequence, tContrast] = GetSequenceContrast(experiment,rule,event_list,tResponse)

% Created 8/28/14 by DJ.

if nargin<4
    tResponse = 0:100:1000;
end

% rule = 'T0vD0';
if strcmp(experiment,'sq')
    tEvents = 0:1000:3750;
else
    tEvents = 0:500:2500;
end


if strcmp(experiment,'sf3')
    switch rule
        case 'T1vT0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pT_{0/3}','pD^*_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pD^*_{-/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdTdd-ddTdd'};
        case 'T2vT1'
            sequence{1} = {'pD_{0/3}','pT_{0/3}','pT_{1/3}','pD_{2/3}','pD^*_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pD_{+/3}','pD_{+/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'ttTdd-dtTdd'};
        case 'T2vT0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pT_{0/3}','pD^*_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pD_{+/3}','pD_{+/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'ttTdd-ddTdd'};            
            
        case 'T0vD0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pD_{0/3}','pD_{0/3}','pT_{0/3}','pD^*_{-/3}','pD_{-/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'ddTdd-dDddd'};
        case 'D1vD0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pD_{1/3}','pD_{1/3}','pD^*_{-/3}','pD_{-/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdDdd-dDddd'};
        case 'T1vD0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pD_{1/3}','pT_{1/3}','pD_{2/3}','pD^*_{-/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdTdd-dDddd'};
        case 'D2vD0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pT_{1/3}','pD_{2/3}','pD_{2/3}','pD^*_{-/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(4);
            legendstr = {'ttDdd-dDddd'}; 
            
        case 'D2vD1'
            sequence{1} = {'pT_{0/3}','pD_{1/3}','pD_{1/3}','pD^*_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pT_{1/3}','pD_{2/3}','pD_{2/3}','pD^*_{-/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(4);
            legendstr = {'ttdDd-tdDdd'};     
        case {'T2vD0','T*vD0'}
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pD_{+/3}','pD_{+/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'ttTdd-dDddd'}; 
            
        case 'T+vD0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pD_{+/3}','pT_{+/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(5);
            legendstr = {'tttdT-dDddd'};  
        case 'D+vD0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pD_{+/3}','pD_{+/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(5);
            legendstr = {'tttdD-dDddd'};    
        case 'T*vD*'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pT_{0/3}','pT_{1/3}','pT^*_{2/3}','pD_{+/3}','pD_{+/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'ttTdd-ddDdd'}; 
        case 'D*vD0'
            sequence{1} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            sequence{2} = {'pD_{0/3}','pD_{0/3}','pD^*_{-/3}','pD_{-/3}','pD_{-/3}','sf3-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'ttTdd-ddDdd'};             
            
        otherwise
            sequence = {{},{}};
    end
else


    switch rule
        case 'TDvDD'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','NULL','sf-Circle'};
            sequence{2} = {'pT_{0/2}','pD_{1/2}','pD_{1/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};        
            tContrast{1} = tResponse + tEvents(1);
            tContrast{2} = tResponse + tEvents(1);
            legendstr = {'Tdddd-Ddddd'};
        case 'T1vT0'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pT_{0/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};
            sequence{2} = {'pT_{0/2}','pD_{1/2}','pT^*_{1/2}','pD_{+/2}','pD_{+/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdTdd-ddTdd'};
        case 'T1vD0'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','NULL','sf-Circle'};
            sequence{2} = {'pT_{0/2}','pD_{1/2}','pT^*_{1/2}','pD_{+/2}','pD_{+/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdTdd-ddDdd'};
        case 'T0sn4vT0sn2'
            sequence{1} = {'pD_{0/2}','pT_{0/2}','pD_{1/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};
            sequence{2} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pT_{0/2}','pD^*_{-/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(4);
            legendstr = {'dddTd-dTddd'};
        case 'T1sn4vT1sn2'
            sequence{1} = {'pT_{0/2}','pT^*_{1/2}','pD_{+/2}','pD_{+/2}','pD_{+/2}','sf-Circle'};
            sequence{2} = {'pD_{0/2}','pD_{0/2}','pT_{0/2}','pT^*_{1/2}','pD_{+/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(2);
            tContrast{2} = tResponse + tEvents(4);
            legendstr = {'ddtTd-tTddd'};
        case 'D1vD0'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','NULL','sf-Circle'};
            sequence{2} = {'pT_{0/2}','pD_{1/2}','pD_{1/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdDdd-ddDdd'};
        case 'T0vD0'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','NULL','sf-Circle'};
            sequence{2} = {'pD_{0/2}','pD_{0/2}','pT_{0/2}','pD_{1/2}','pD^*_{-/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'ddTdd-ddDdd'};
            
        case 'T+vD0'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','NULL','sf-Circle'};
            sequence{2} = {'pT_{0/2}','pT^*_{1/2}','pD_{+/2}','pT_{+/2}','pD_{+/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(4);
            legendstr = {'ttdTd-ddDdd'};   
        case 'D+vD0'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','NULL','sf-Circle'};
            sequence{2} = {'pT_{0/2}','pT^*_{1/2}','pD_{+/2}','pD_{+/2}','pD_{+/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(4);
            legendstr = {'ttdDd-ddDdd'};  
        
        case 'T*vD*'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','pD_{-/2}','sf-Circle'};
            sequence{2} = {'pT_{0/2}','pD_{1/2}','pT^*_{1/2}','pD_{+/2}','pD_{+/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(4);
            tContrast{2} = tResponse + tEvents(3);
            legendstr = {'tdTdd-dddDd'};
            
        case 'D*vD0'
            sequence{1} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','pD_{-/2}','sf-Circle'};            
            sequence{2} = {'pD_{0/2}','pD_{0/2}','pD_{0/2}','pD^*_{-/2}','pD_{-/2}','sf-Circle'};
            tContrast{1} = tResponse + tEvents(3);
            tContrast{2} = tResponse + tEvents(4);
            legendstr = {'dDddd-dddDd'};  
        otherwise
            sequence = {{},{}};
    end
end

if strcmp(experiment,'sq')
    for i=1:2
        for j=1:numel(sequence{i})
            if sequence{i}{j}(1) == 'p'
                sequence{i}{j}(1) = 'a';
            end
            if strcmp(sequence{i}{j}(1:3),'sf-') % for sf-circle events
                sequence{i}{j}(1:3) = '';
            end
        end
    end
end

C_minus = BuildSequenceContrast(event_list,tResponse,sequence{1},tEvents,tContrast{1});
C_plus = BuildSequenceContrast(event_list,tResponse,sequence{2},tEvents,tContrast{2});

contrasts = C_plus-C_minus;
% Plot contrast matrix
figure(55);
imagesc(tResponse,1:size(contrasts,1),contrasts);
xlabel('tResponse')
ylabel('Event')
T = length(tResponse);
set(gca,'ytick',(1:T:size(contrasts,1)),'yticklabel',event_list);