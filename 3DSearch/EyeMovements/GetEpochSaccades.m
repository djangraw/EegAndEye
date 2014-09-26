function saccade_times = GetEpochSaccades(EEG, saccade_codes, t_discrim)

% Finds the times at which saccades were made in each epoch.
%
% saccade_times = GetEpochSaccades(EEG, saccade_codes, t_discrim)
%
% INPUTS:
% -EEG is an epoched EEGLAB data structure from a 3DSearch experiment.
% -saccade_codes is a vector of event numbers whose times you want to find 
%  in each epoch, or a cell array of strings indicating the event types 
%  (see GetNumbers for the possible constants and their names).  
%  [Default is 'SACCADE_END'.]
% -t_discrim is a vector containing the time range (in ms) relative to
%  t=0 in which saccades will be considered potential discriminating
%  saccades. (default: [-inf inf])
%
% OUTPUTS:
% -saccade_times is a cell array with length equal to the number of epochs
%  in the EEG struct.  Each cell contains a vector of the times, in ms, at
%  which saccades ended during that epoch.
%
% Created 5/12/11 by DJ.
% Updated 6/2/11 by DJ - switched from SACCADE_TO to SACCADE_END event code
% Updated 6/10/11 by DJ - added optional saccade_codes input
% Updated 8/2/11 by DJ - added t_discrim input

if nargin<2
    saccade_codes = 'SACCADE_END';
end
if nargin<3
    t_discrim = [-inf inf]; % epoch time range (min, max in ms) in which saccades could be used in discrimination.
%     t_discrim = [0 1000]; 
end

% convert string or numeric inputs to corresponding codes
if ischar(saccade_codes) % e.g. 'SACCADE_START'
    GetNumbers;
    saccade_codes = {num2str(Numbers.(saccade_codes))};
elseif isnumeric(saccade_codes) % e.g. [2000 2001]
    tmp = saccade_codes;
    saccade_codes = cell(0);
    for i=1:numel(tmp)
        saccade_codes{i} = num2str(tmp(i));
    end
elseif iscell(saccade_codes) % e.g. {'SACCADE_START', 'SACCADE_END'}
    if ischar(saccade_codes{1})
        GetNumbers;
        for i=1:numel(saccade_codes)
            saccade_codes{i} = num2str(Numbers.(saccade_codes{i}));
        end
    end
end

% set up
nTrials = numel(EEG.epoch);
saccade_times = cell(nTrials,1);

% find saccade times
for i=1:nTrials
    % Find times of all saccades in epoch
    isSaccadeEvent = ismember(EEG.epoch(i).eventtype,saccade_codes); % find the events whose labels match the saccade code
    allSaccadeTimes = [EEG.epoch(i).eventlatency{isSaccadeEvent}]; % get the corresponding times for these events

    % Find times of saccades in acceptable discrimination window
    okSaccadeTimes = allSaccadeTimes(allSaccadeTimes>=min(t_discrim) & allSaccadeTimes<=max(t_discrim)); 
    % Make sure you have at least one
    if isempty(okSaccadeTimes) % if not...
        okSaccadeTimes = allSaccadeTimes(find(allSaccadeTimes>=min(t_discrim),1,'first')); % Get first saccade after min(t_discrim)
    end
    
    % Pass result to output variable
    saccade_times{i} = okSaccadeTimes;
    
end