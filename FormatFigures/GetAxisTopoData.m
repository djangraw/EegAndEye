function [data,chanx,chany] = GetAxisTopoData(hAxis)

% Create 5/8/13 by DJ.

kids = get(hAxis,'Children');

% Get channel locations
types = get(kids,'Type');
linekids = kids(strcmp('line',types));
markers = get(linekids,'Marker');
channelkid = linekids(strcmp('.',markers));
chanx = get(channelkid,'XData');
chany = get(channelkid,'YData');

%% Get topo data
surfacekid = kids(strcmp('surface',types));
surfacex = get(surfacekid,'XData');
surfacey = get(surfacekid,'YData');
surfacec = get(surfacekid,'CData');

ok = ~isnan(surfacec);
% data = griddata(surfacex(ok),surfacey(ok),surfacec(ok),chanx,chany,'invdist');
data = griddata(surfacex(ok),surfacey(ok),surfacec(ok),chanx,chany,'linear');
% data = TriScatteredInterp(surfacex(ok),surfacey(ok),surfacec(ok),chanx,chany,'invdist');