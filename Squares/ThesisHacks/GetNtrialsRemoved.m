% function GetNtrialsRemoved
% Created 4/23/14 by DJ.

clear nRemoved*

R = R_sf3_sqnum;
trial_rej_rules = {'early_button' 'late_button' 'wrong_button'};
subjects = 1:12;%[1:10 12:13];
prefix = 'sf3';

% R = R_sq_sqnum;
% trial_rej_rules = {'skipped_ends'  'skip'  'backward' 'early_button' 'late_button' 'wrong_button'};
% subjects = [9:11, 13:15 17:27];
% prefix = 'sq';

for i=1:numel(R)
    y = loadBehaviorData(subjects(i),[],prefix);
    EEG = R.EEG;
    EEG.etc.rejectepoch(:) = false;
    EEG = RejectEegData(EEG,y,trial_rej_rules);
    
    nRemoved(i) = sum(R(i).EEG.etc.rejectepoch);
    nRemovedForBehavior(i) = sum(EEG.etc.rejectepoch);
    nRemovedForVoltage(i) = nRemoved(i)-nRemovedForBehavior(i);
end
fprintf('%g for behavior, %g for voltage\n', mean(nRemovedForBehavior),mean(nRemovedForVoltage));

