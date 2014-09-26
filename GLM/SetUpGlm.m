function [EEG, regressor_events, artifact_events, nuisance_events] = SetUpGlm(subject,EEG,method,offset,rejectionrule,extent_ms,artifact_extent_ms)

% [EEG, regressor_events, artifact_events, nuisance_events] = SetUpGlm(subject,EEG,method,offset,rejectionrule,extent_ms,artifact_extent_ms)
%
% INPUTS:
% - subject is a scalar indicating the subject number whose behavioral data
%  you would like to use.
% - EEG is the subject's corresponding EEG dataset (all subjects combined).
% - method is a string indicating which kind of regressors you want to use.
%  (see code for options.)
% - offset is a scalar indicating the offset, in ms, of events.  This is a 
%  manual correction for the misalignment sometimes seen in the 
%  Eyelink/Sensorium pairing).  [default = 0].
% - rejectionrule is a string or cell array of strings indicating which
%  RejectBehavior rules should be used to reject trials.
% - extent_ms and artifact_extent_ms are scalars or 2-element vectors
% indicating how many samples will be influenced by each event or artifact.
%
% OUTPUTS:
% - EEG is an EEGLAB data struct, with the extra events added.
% - regressor_events is a cell array of strings indicating the regressor
%  event names.
% - artifact_events is a cell array of strings indicating the artifact
%  event names.
% - nuisance_events is a cell array of strings indicating the nuisance
% regressor event names.
% 
% Created 12/23/11 by DJ.
% Updated 1/4/11 by DJ - added method input.
% Updated 1/10-12/12 by DJ - added more options for method
% Updated 2/1/12 by DJ - added nuisance_events
% Updated 2/8/12 by DJ - added -Start option
% Updated 3/22/12 by DJ - added Saccade-SqNum option
% Updated 4/9/12 by DJ - added v1pt5 options
% Updated 4/20/12 by DJ - added offset input
% Updated 6/11/12 by DJ - added SqNum-v1pt5, Type-v1pt5 options
% Updated 6/15/12 by DJ - added Saccade-LogOdds-v1pt5 option
% Updated 7/20/12 by DJ - added eventweights
% Updated 7/30/12 by DJ - sped up by adding events all at once
% Updated 8/8/12 by DJ - switched to ureventweights, added epoching code
% Updated 3/22/13 by DJ - added artifact_extent_ms input

% Handle inputs
if nargin<3 || isempty(method)
    method = '9reg';
end
if nargin<4 || isempty(offset)
    offset = 0;
end
if nargin<5
    rejectionrule = {'backward' 'skipped_square' 'skipped_ends'};
end
if nargin<6 || isempty(extent_ms)
    extent_ms = 500;
end
if nargin<7 || isempty(artifact_extent_ms)
    artifact_extent_ms = 500;
end
if ~iscell(rejectionrule)
    rejectionrule = {rejectionrule};
end
vThreshold = 75; % max voltage allowed near a regressor/nuisance event, in uV

% Get event types
switch method
    case 'Saccade-DiffLogOdds-v1pt5'
        regressor_events = {'DiffLogOddsSquare','DiffLogOddsCompl','DiffLogOddsIncompl','DiffLogOddsIrrel','Circle-D','Circle-T'};
        artifact_events = {'Cross','Errant','BlinkStart','BlinkEnd','Button'};
        nuisance_events = {'Saccade'};
    case 'Saccade-LogOdds-v1pt5'
        regressor_events = {'LogOddsDist','LogOddsTarg','LogOddsIrrel','Circle-D','Circle-T'};
        artifact_events = {'Cross','Errant','BlinkStart','BlinkEnd','Button'};
        nuisance_events = {'Saccade'};
    case 'SqNum-v1pt5'
        regressor_events = {'SqNum1', 'SqNum2','SqNum3','SqNum4','SqNum5','Circle-D','Circle-T'};
        artifact_events = {'Cross','Errant','BlinkStart','BlinkEnd','Button'};
        nuisance_events = {};
    case 'Type-v1pt5'
        regressor_events = {'Dist-0T','Dist-1T','Integ','Compl','Irrel','Circle-D','Circle-T'};
        artifact_events = {'Cross','Errant','BlinkStart','BlinkEnd','Button'};
        nuisance_events = {};
    case 'SqNum-Type-v1pt5'
        regressor_events = {'Dist-0T','Dist-1T','Integ','Compl','Irrel'};
        artifact_events = {'Cross','Errant','BlinkStart','BlinkEnd','Button'};
        nuisance_events = {'SqNum1', 'SqNum2','SqNum3','SqNum4','SqNum5','Circle-D','Circle-T'};
    case 'Type-SqNum-v1pt5'
        regressor_events = {'SqNum1', 'SqNum2','SqNum3','SqNum4','SqNum5'};
        artifact_events = {'Cross','Errant','BlinkStart','BlinkEnd','Button'};
        nuisance_events = {'Dist-0T','Dist-1T','Integ','Compl','Irrel','Circle-D','Circle-T'};
    case 'Saccade-LastSq-v1pt5'
        regressor_events = {'Sq1-4','SqNum5','Circle-D','Circle-T'};
        artifact_events = {'Cross','Errant','BlinkStart','BlinkEnd','Button'};
        nuisance_events = {'Saccade'};
    case 'Saccade-Type-v1pt5'
        regressor_events = {'Dist-0T','Dist-1T','Integ','Compl','Irrel','Circle-D','Circle-T'};
        artifact_events = {'Cross','Errant','BlinkStart','BlinkEnd','Button'};
        nuisance_events = {'Saccade'};
    case 'Saccade-SqNum-v1pt5'
        regressor_events = {'SqNum1', 'SqNum2','SqNum3','SqNum4','SqNum5','Circle-D','Circle-T'};
        artifact_events = {'Cross','Errant','BlinkStart','BlinkEnd','Button'};
        nuisance_events = {'Saccade'};
    case 'Saccade-SqNum' % square number
        regressor_events = {'SqNum1', 'SqNum2','SqNum3','SqNum4','SqNum5','Circle-D','Circle-T','Button-D','Button-T'};
        artifact_events = {'Cross','Errant', 'BlinkStart','BlinkEnd'};
        nuisance_events = {'Saccade'};
    case 'Saccade-9reg' % 9 reg plus Saccade
        regressor_events = {'Dist-0T','Dist-1T','Integ','Compl','Irrel','Circle-D','Circle-T','Button-D','Button-T'};
        artifact_events = {'Cross','Errant', 'BlinkStart','BlinkEnd'};
        nuisance_events = {'Saccade'};
    case 'Saccade-9reg-Start' % 9 reg plus Saccade (2/8/12)
        regressor_events = {'SacStartDist-0T','SacStartDist-1T','SacStartInteg','SacStartCompl','SacStartIrrel','SacStartCircle-D','SacStartCircle-T','Button-D','Button-T'};
        artifact_events = {'SacStartCross','SacStartErrant', 'BlinkStart','BlinkEnd'};
        nuisance_events = {'StartSaccade'};    
    case '9reg' % 9 regressors (1/4/12)
        regressor_events = {'Dist-0T','Dist-1T','Integ','Compl','Irrel','Circle-D','Circle-T','Button-D','Button-T'};
        artifact_events = {'Cross','Errant', 'BlinkStart','BlinkEnd'};
        nuisance_events = {};
    case '7reg' % 7 regressors (12/23/11)
        regressor_events = {'Dist-0T','Dist-1T','Dist-2T','Integ','Compl','Button-D','Button-T'};
        artifact_events = {'BlinkStart','BlinkEnd'};
        nuisance_events = {};
    case '10con' % 10 contrasts (1/10/12) - (Non-orthogonal, doesn't work well)
        regressor_events = {'Prep','Antic','Dist','Targ','Integ','Compl','Square','Saccade','Circle','Button'};
        artifact_events = {'Cross','Errant', 'BlinkStart','BlinkEnd'};
        nuisance_events = {};
    case '6con' % 6 contrasts (1/11/12) - (Non-orthogonal, doesn't work well)
        regressor_events = {'Prep','Antic','Dist','Targ','Saccade','Button'};
        artifact_events = {'Cross','Errant', 'BlinkStart','BlinkEnd'};
        nuisance_events = {};
    case '4con' % 4 contrasts (1/11/12) - (Non-orthogonal, doesn't work well)
        regressor_events = {'Prep','Antic','Dist','Targ'};
        artifact_events = {'Cross','Errant', 'BlinkStart','BlinkEnd'};    
        nuisance_events = {};
    otherwise
        error('Method not recognized!')
end

disp('Getting event times...');
% load behavioral data
y = loadBehaviorData(subject);
% find events
all_events = [regressor_events artifact_events nuisance_events]; 
nEvents = numel(all_events);
nSesssions = numel(y);
[times, codes, weights] = deal(cell(nEvents,nSesssions)); % events by sessions
for i=1:nEvents
    if ~any(strcmp(all_events{i},{EEG.event.type})) % if this event hasn't been added yet...
        [times(i,:),codes(i,:),weights(i,:)] = UseEventRule(y,all_events{i});        
        for j=1:nSesssions
            times{i,j} = times{i,j} + offset; % add offset to times
        end        
    end
end
% append events for each session
[all_times, all_codes] = deal(cell(1,nSesssions));
for j=1:nSesssions
    all_times{j} = cat(1,times{:,j});
    all_codes{j} = cat(1,codes{:,j});    
end

% Add events
disp('Adding events to data struct...');
EEG = AddEeglabEvents_MultiSession(EEG,y,all_times,all_codes);

% Once all the events are added, add the weights
if ~isfield(EEG.etc,'ureventweights')
    EEG.etc.ureventweights = ones(size(EEG.event));
end
for i=1:nEvents
    if ~isempty(weights{i})
        EEG.etc.ureventweights(strcmp(all_events{i},{EEG.urevent(:).type})) = cat(1,weights{i,:});
    end
end

% Epoch data
EEG = pop_epoch(EEG,{'TrialStart-T', 'TrialStart-D'},[-1.5 4.5]);
% Interpolate noisy electrodes and reject noisy trials
EEG = EnforceVoltageThreshold(EEG,vThreshold,[nuisance_events regressor_events], artifact_events, extent_ms, artifact_extent_ms);

% Reject bad trials
disp(['Rejecting trials according to rules: ' sprintf('%s, ', rejectionrule{:}) '...']);
EEG = RejectEegData(EEG,y,rejectionrule); % Rename as 'REJECTED' any events in a bad range

disp('Done!---Now run the folowing line of code:')
disp('[responseFns, tResponse, NewEEG] = RunEegGlm(EEG,regressor_events,extent_ms,artifact_events,artifact_extent_ms);')