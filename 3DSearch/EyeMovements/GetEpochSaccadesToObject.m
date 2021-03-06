function [saccade_start_times, saccade_end_times] = GetEpochSaccadesToObject(BigSaccade,EEG,dist_threshold,t_discrim)

% Finds the times at which saccades were made to the object or away from it
% in each epoch.
%
% [saccade_start_times, saccade_end_times] = GetEpochSaccadesToObject(BigSaccade,EEG,dist_threshold,t_discrim) 
%
% Like GetEpochSaccades, but (1) works off of BigSaccade instead of EEG
% event info, and (2) allows the use of a distance threshold.
%
% INPUTS:
% -BigSaccade is a bigsaccade struct created using GetBigSaccadeStruct
% (called with the same EEG input as this function's EEG input).
% -EEG is an epoched EEGLAB data structure from a 3DSearch experiment.
% -dist_threshold is a scalar indicating the max distance from an object
% that will still count as a saccade to the object.
% -t_discrim is a vector containing the time range (in ms) relative to
%  t=0 in which saccades will be considered potential discriminating
%  saccades. (default: [-inf inf])
%
% OUTPUTS:
% -saccade_start_times is a cell array with length equal to the number of epochs
%  in the EEG struct.  Each cell contains a vector of the times, in ms, at
%  which saccades to or from the object started during that epoch.
% -saccade_end_times is a cell array with length equal to the number of epochs
%  in the EEG struct.  Each cell contains a vector of the times, in ms, at
%  which saccades to or from the object ended during that epoch.
%
%
% Created 8/16/11 by DJ
 
% Set up
if nargin<4
    t_discrim = [-inf inf]; % epoch time range (min, max in ms) in which saccades could be used in discrimination.
%     t_discrim = [0 1000]; 
end


setnumber = find(strcmp(EEG.setname,BigSaccade.EEGsetnames));

if isempty(setnumber)
    error('This data structure was not found in the BigSaccade struct!');
end

saccade_start_times = cell(EEG.trials,1);
saccade_end_times = cell(EEG.trials,1);
isToObject = (BigSaccade.saccade_end_disttoobj<dist_threshold);
isFromObject = (BigSaccade.saccade_start_disttoobj<dist_threshold & BigSaccade.saccade_end_disttoobj>dist_threshold);

% Find them!
for i=1:EEG.trials
    % get indices of BigSaccade that are in this epoch 
    isInEpoch = (BigSaccade.EEGset==setnumber & BigSaccade.EEGepoch==i);  
    
    % include saccades to object and saccades away from object
    candidateSaccades = find(isInEpoch & (isToObject | isFromObject));
    candidateSaccadeTimes = BigSaccade.saccade_end_latency(candidateSaccades); 
    
    % Find times of saccades in acceptable discrimination window
    isOkCandidate = (candidateSaccadeTimes>=min(t_discrim) & candidateSaccadeTimes<=max(t_discrim)); 
    % Make sure you have at least one
    if sum(isOkCandidate)==0 % if not...
        isOkCandidate = find(candidateSaccadeTimes>=min(t_discrim),1,'first'); % Get first saccade after min(t_discrim)
    end
    
    % Pass result to output variable
    saccade_end_times{i} = candidateSaccadeTimes(isOkCandidate)'; % make row vector of saccade end times
    saccade_start_times{i} = BigSaccade.saccade_start_latency(candidateSaccades(isOkCandidate))'; % make row vector of saccade start times
end