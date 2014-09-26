function y = AddSyncField(y0)

% y = AddSyncField(y0)
%
% Add a sync field to a 3DS data struct to let us use the
% AddEeglabEvents_MultiSession function.
%
% INPUTS:
% -y0 is a vector of 3DS behavior structs.
% 
% OUTPUTS:
% -y is the same as the input, but with sync fields added to each struct.
% 
% Created 6/4/13

y = y0; % start with original
if isfield(y0(1),'sync')
    disp('sync field is already present!')    
    return
else
    
    % add sync field
    for i=1:numel(y)
        y(i).sync.eyelink = y0(i).eyelink.sync_events(:,1);
        y(i).sync.eeg = y0(i).eeg.sync_events(:,1);
        y(i).sync.events = y0(i).eyelink.sync_events(:,2);
    end
end