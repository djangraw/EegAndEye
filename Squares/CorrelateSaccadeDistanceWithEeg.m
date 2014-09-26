function CorrelateSaccadeDistanceWithEeg(y,EEG,eventtype,trial_rej_rules)

% Created 8/2/12 by DJ.

% Get EEG epochs
[~, times, epochs] = GetSquaresErps(EEG,{eventtype},[-498 500]);

% Get saccade distance

% Find relevant saccades
    
s = combine_structs([y.saccade],[],'columns');
[~,~,~,iGood] = UseEventRule(y,eventtype);
for j=1:numel(y)
    [~,isBad] = RejectBehavior(y(j),trial_rej_rules);
    iEvents{j} = setdiff(iGood{j},find(isBad));
end
iOk = iEvents{1};
offset = 0;
for j=2:numel(iEvents)
    offset = offset + numel(y(j-1).saccade.end_time);
    iOk = [iOk; iEvents{j}+offset];
end

% Get saccade info
distance = sqrt((s.end_position(iOk,1)-s.start_position(iOk,1)).^2 + ...
    (s.end_position(iOk,2)-s.start_position(iOk,2)).^2);

if numel(distance)~=size(epochs{1},3)
    error('Number of epochs don''t match!')
end

meanAcrossElectrodes = squeeze(mean(epochs{1},1));
[~,iMax] = max(abs(meanAcrossElectrodes));

scatter(distance,meanAcrossElectrodes(iMax));

