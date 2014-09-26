function comp = GetBlinkComponent(EEG,offset_ms,blinktypes)

% comp = GetBlinkComponent(EEG,offset_ms,blinkevents)
%
% INPUTS:
% -EEG is an EEGLAB dataset with events of type 'BS' and 'BE' added (to
% signify the start and end of each blink).
% -offset_ms is a scalar indicating the offset, in ms, from the stated
% event times when the event actually occurred.
% -blinktypes is a 2-element cell array of strings indicating the event
% markers that signify the start and end of blinks in the EEG struct.
%
% OUTPUTS:
% -comp is a Dx1 vector indicating the component showing maximum power
% during blink periods.
%
% Created 6/4/13 by DJ.
% Updated 3/11/14 by DJ - added check to remove unfinished blinks, commments.
% Updated 9/26/14 by DJ - added blinktypes input.

if nargin<2 || isempty(offset_ms)
    offset_ms = 0;
end
if nargin<3 || isempty(blinktypes)
    blinktypes = {'BS','BE'};
end

% Set options
doPlot = false;

% Get times of blink start and end
eventType = {EEG.event.type};
eventLatency = [EEG.event.latency];
iEventsBS = find(strcmp(blinktypes{1},eventType));
iEventsBE = find(strcmp(blinktypes{2},eventType));
iEventsBound = find(strcmp('boundary',eventType));

% Remove unfinished blinks
for i=1:numel(iEventsBound)
    if iEventsBS(find(iEventsBS<iEventsBound(i),1,'last')) > iEventsBE(find(iEventsBE<iEventsBound(i),1,'last'))
        iEventsBS(find(iEventsBS<iEventsBound(i),1,'last')) = [];
        fprintf('Ignoring blink start at end of session %d\n',i)
    end
end
if iEventsBS(end) > iEventsBE(end)
    iEventsBS(end) = [];
end
fprintf('%d blink starts, %d blink ends\n',length(iEventsBS),length(iEventsBE));

% Set up artifact checks
energy = mean(EEG.data.^2,1); % mean energy across all electrodes
ENERGY_MAX = 1e4; % If the energy is greater than this, we're in an unusual artifact, not a regular blink.
skipped = 0;
% Make vector of wheter blinks were in progress
isBlink = false(1,EEG.pnts);
offset_samples = round(offset_ms/1000*EEG.srate);
for i=1:numel(iEventsBS)
    iBlink = max(1, round(eventLatency(iEventsBS(i))+offset_samples)) : min(EEG.pnts, round(eventLatency(iEventsBE(i))+offset_samples));
    if max(energy(iBlink)) < ENERGY_MAX % if it's not a weird artifact
        isBlink(iBlink) = true;
    else
        skipped = skipped + 1;
    end
end
fprintf('Skipped %d/%d blinks due to artifacts.\n',skipped,numel(iEventsBS));

% Extract data from these times
data_blink = EEG.data(:,isBlink);

% Find max power component
comp = Maximumpower(data_blink);

% plot component on scalp
if doPlot
    cla;
    topoplot(comp,EEG.chanlocs);
    colorbar;
    title(sprintf('Blink compoent for %s',EEG.setname));
end

% Function from eyesubtract.m
function Max_eigvec=Maximumpower(EEGdata)

[Channel,~,Epoch]=size(EEGdata);

for e=1:Epoch
    for i=1:Channel
        EEGdata(i,:,e)=EEGdata(i,:,e)-mean(EEGdata(i,:,e));
    end;
end;


% Maximum Power method

a=EEGdata(:,:);

pow_a=a*a';

[vb,tmp] = eig(pow_a);
Max_eigvec=vb(:,end)./norm(vb(:,end));


return