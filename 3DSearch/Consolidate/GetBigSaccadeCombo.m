function Acombo = GetBigSaccadeCombo(EEGdis,EEGtarg,x,pt)

% Get a saccade info struct that combines both target and distractor info.
%
% Acombo = GetBigSaccadeCombo(EEGdis,EEGtarg,x,pt)
%
% INPUTS:
% - EEGdis and EEGtarg are eeglab datasets containing all distractor epochs 
% and target epochs, respectively.
% - x is a vector of 3DSearch behavior structs as imported using
% Import_3DS_Data_v3
% - pt is a matrix of posterior probabilities of this saccade being the
% "aha" saccade, as calculated by our jittered logistic regression
% algorithm and loaded using GetFinalPosteriors.
% 
% OUPTUTS:
% - Acombo is a struct containing info about every potential "aha" saccade 
% in every trial.  Start and end time and position, eeg epoch number,
% posterior probability, etc. are included.  Each field is a column vector
% in which each row represents a saccade.
%
% Created 7/25/11 by DJ.
% Updated 7/27/11 by DJ - added use of first saccade after time range if 
% none are found in the range

% Get info
Adis = GetBigSaccadeStruct(EEGdis,x);
Atarg = GetBigSaccadeStruct(EEGtarg,x);

% Crop
% Adis_cropped = crop_struct(Adis,~isnan(Adis.EEGepoch));
% Atarg_cropped = crop_struct(Atarg,~isnan(Atarg.EEGepoch));

% Add posterior field
Adis.posterior = nan(size(Adis.EEGepoch));
Atarg.posterior = nan(size(Atarg.EEGepoch));

for i=1:EEGdis.trials
    % find saccades and posteriors in range
    okPosteriors = pt(i,pt(i,:)>0);
    okSaccades = find(Adis.EEGepoch==i & Adis.saccade_end_latency>0 & Adis.saccade_end_latency<1000);
    % as in the algorithm, if no saccades are within the time range we use the first one after the range.
    if isempty(okSaccades)
        okSaccades = find(Adis.EEGepoch==i & Adis.saccade_end_latency>0,1);
    end
    % consistency check
    if numel(okPosteriors) ~= numel(okSaccades)
        warning('GetBigSaccadeCombo:epochMismatch','Distractor epoch %d doesn''t match\n',i)
    else
        % Add posterior to struct
        Adis.posterior(okSaccades,:) = okPosteriors';
    end
end
for i=1:EEGtarg.trials
    % find saccades and posteriors in range
    okPosteriors = pt(EEGdis.trials+i,pt(EEGdis.trials+i,:)>0);
    okSaccades = find(Atarg.EEGepoch==i & Atarg.saccade_end_latency>0 & Atarg.saccade_end_latency<1000);
     % as in the algorithm, if no saccades are within the time range we use the first one after the range.
    if isempty(okSaccades)
        okSaccades = find(Atarg.EEGepoch==i & Atarg.saccade_end_latency>0,1);
    end
    % consistency check
    if numel(okPosteriors) ~= numel(okSaccades)
        warning('GetBigSaccadeCombo:epochMismatch','Target epoch %d doesn''t match\n',i)
    else
        % Add posterior to struct
        Atarg.posterior(okSaccades,:) = okPosteriors';
    end
end

% crop
Adis_cropped = crop_struct(Adis,~isnan(Adis.posterior));
Atarg_cropped = crop_struct(Atarg,~isnan(Atarg.posterior));

% combine
Acombo = combine_structs(Adis_cropped,Atarg_cropped,'columns');

% add isTargetEpoch field
Acombo.isTargetEpoch = [zeros(size(Adis_cropped.EEGepoch)); ones(size(Atarg_cropped.EEGepoch))];
