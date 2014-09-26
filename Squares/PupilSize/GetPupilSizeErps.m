function ps_erp = GetPupilSizeErps(prefix,subject,all_events,tEpoch,tBaseline)

if numel(tEpoch)==2 % convert from limits to vector
    tEpoch = tEpoch(1):tEpoch(2);
end
if numel(tBaseline)==2 % convert from limits to vector
    tBaseline = tBaseline(1):tBaseline(2);
end

% Load behavioral data
load(sprintf('%s-%d-all.mat',prefix,subject));
nSessions = numel(y);
% Load pupil size
ps = cell(1,nSessions);
for j=1:nSessions
    foo = load(sprintf('%s-%d-%d-eyepos.mat',prefix,subject,y(j).session));
    % Interpolate Blinks    
    tPs = (1:length(foo.pupilsize)) + y(j).recording_start_time - 1;
    ps{j} = InterpolateBlinks(foo.pupilsize,tPs,y(j));
end
% % Get event times/IDs
% events = cell(1,nSessions);
% for j=1:nSessions
%     events{j} = [y(j).saccade.end_time, y(j).saccade.class];
%     events{j} = events{j}(~isnan(y(j).saccade.class),:);
% end

% Get different event types
nEvents = numel(all_events);
[times, codes, weights] = deal(cell(nEvents,nSessions)); % events by sessions
for i=1:nEvents
    [times(i,:),codes(i,:),weights(i,:)] = UseEventRule(y,all_events{i});        
    for j=1:nSessions
        times{i,j} = times{i,j} - y(i).recording_start_time + 1; % add offset to times
    end        
end 



% Get avg activity
nT = length(tEpoch);
ps_erp = zeros(nEvents,nT);
ps_epoch = cell(nEvents,1);
ps_epoch_session = cell(nEvents,nSessions);
for i=1:nEvents
    for j=1:nSessions
        % Get epochs for this session & event type
        ps_epoch_session{i,j} = nan(length(times{i,j}),nT);
        for k=1:length(times{i,j})
            try
                ps_epoch_session{i,j}(k,:) = ps{j}(times{i,j}(k) + tEpoch);
            catch % if we exceed the bounds
                iOk = (times{i,j}(k) + tEpoch) > 0 & (times{i,j}(k) + tEpoch) < length(ps{j});
                ps_epoch_session{i,j}(k,iOk) = ps{j}(times{i,j}(k) + tEpoch(iOk));
            end
        end                
    end
    % Take average across sessions
    ps_epoch{i} = cat(1,ps_epoch_session{i,:});
    ps_erp(i,:) = nanmean(ps_epoch{i},1);
    % Subtract baseline
    ps_erp(i,:) = ps_erp(i,:) - nanmean(ps_erp(i,tBaseline+1));
end
