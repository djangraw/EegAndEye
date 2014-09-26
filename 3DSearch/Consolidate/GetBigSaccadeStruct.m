function A = GetBigSaccadeStruct(EEG,x,sessionNumber)

% Create a struct with all the saccade info in it.
%
% A = GetBigSaccadeStruct(EEG,x,sessionNumber)
%
% INPUTS: 
% - EEG is an eeglab data struct from a 3DSearch or Squares experiment or 
% an m-element vector of such structs.
% - x is a 3DSearch or Squares behavior data struct imported using 
% Import_3DS_Data_v3/import_squares_data or an n-element vector of such 
% structs.
% - sessionNumber is a scalar or n-element vector of the session number (if
% EEG contains multiple sessions). [1=first or only session]
%
% OUTPUTS:
% - A is a struct with fields relating to the properties of each saccade.
% For example, A.saccade_start_time(i) is the start time of saccade #i.
%
% Created 7/22/11 by DJ.
% Updated 7/25/11 by DJ - resolved session ambiguities, added multi-session
% capabilities
% Updated 8/16/11 by DJ - added disttoobj data, multiple EEG file support
% Updated 11/17/11 by DJ - added squares data file support

% Handle defaults
if nargin<3
    sessionNumber = 1:numel(x);
end

% INITIALIZE
A = struct('EEGset',[],'EEGepoch',[],...
    'saccade_start_time',[],'saccade_end_time',[],...
    'saccade_start_latency',[],'saccade_end_latency',[],...
    'saccade_start_x',[],'saccade_start_y',[],...
    'saccade_end_x',[],'saccade_end_y',[],...
    'saccade_duration',[],'saccade_length',[],'saccade_velocity',[],...
    'saccade_start_disttoobj',[],'saccade_end_disttoobj',[],...
    'eyelink_session',[],'eyelink_trial_number',[],...
    'eyelink_trial_time',[],'eyelink_trial_istarget',[],...
    'eyelink_start_latency',[],'eyelink_end_latency',[]);

% Handle multiple x structs
if numel(x)>1
    disp('Getting Big Saccade Struct...'); % this could be a long process, so give status updates
    A = [];
    for i=1:numel(x)
        fprintf('Session %d of %d...\n',i,numel(x));
        Anew = GetBigSaccadeStruct(EEG,x(i),sessionNumber(i));
        A = combine_structs(A,Anew,'columns');        
    end
    A.EEGsetnames = A.EEGsetnames(1,:); % this is the only field that wasn't a column vector
    disp('Success!')
    
else % single x struct

    switch EEG(1).setname(1:2) % what kind of dataset is this?
        case '3D'
            % EXTRACT
            A.saccade_start_time = x.eyelink.saccade_start_times;
            A.saccade_end_time = x.eyelink.saccade_times;
            A.saccade_start_x = x.eyelink.saccade_start_positions(:,1);
            A.saccade_start_y = x.eyelink.saccade_start_positions(:,2);
            A.saccade_end_x = x.eyelink.saccade_positions(:,1);
            A.saccade_end_y = x.eyelink.saccade_positions(:,2);
            % nSaccades = numel(A.saccade_start_time);

            % SORT
            % which saccades are in which epochs? what are their latencies in that epoch?
            [A.EEGset, A.EEGepoch, A.saccade_start_latency, A.saccade_end_latency] = deal(nan(size(A.saccade_start_time)));

            for j=1:numel(EEG)
                sacStartTimes = EyelinkToEegTimes(x.eyelink.saccade_start_times,x)/x.eeg.eventsamplerate;
                [~, saccade_start_latency] = TimeToEpochNumber(EEG(j),sacStartTimes,sessionNumber,1);
                sacEndTimes = EyelinkToEegTimes(x.eyelink.saccade_times,x)/x.eeg.eventsamplerate;
                [EEGepoch, saccade_end_latency] = TimeToEpochNumber(EEG(j),sacEndTimes,sessionNumber,1);

                A.saccade_start_latency(~isnan(EEGepoch)) = saccade_start_latency(~isnan(EEGepoch));
                A.saccade_end_latency(~isnan(EEGepoch)) = saccade_end_latency(~isnan(EEGepoch));
                A.EEGepoch(~isnan(EEGepoch)) = EEGepoch(~isnan(EEGepoch));
                A.EEGset(~isnan(EEGepoch)) = j;
                A.EEGsetnames{j} = EEG(j).setname;
            end

            % CALCULATE
            % duration
            A.saccade_duration = A.saccade_end_time - A.saccade_start_time;
            % distance
            x1 = A.saccade_start_x;
            x2 = A.saccade_end_x;
            y1 = A.saccade_start_y;
            y2 = A.saccade_end_y;
            A.saccade_length = sqrt((x2-x1).^2 + (y2-y1).^2);
            % average velocity
            A.saccade_velocity = A.saccade_length ./ A.saccade_duration;

            % LOCATE
            A.saccade_start_disttoobj = DistToObject(x,[A.saccade_start_x A.saccade_start_y], A.saccade_start_time);
            A.saccade_end_disttoobj = DistToObject(x,[A.saccade_end_x A.saccade_end_y], A.saccade_end_time);

            % TRACK    
            A.eyelink_session = repmat(x.session,size(A.saccade_start_time));
            trial_start_times = x.eyelink.object_events(1:2:end,1);
            trial_end_times = x.eyelink.object_events(2:2:end,1);
            object_numbers = x.eyelink.object_events(1:2:end,2)-500;
            object_istarget = strcmp({x.objects(object_numbers).tag},'TargetObject');

            [A.eyelink_trial_number, A.eyelink_trial_time, A.eyelink_trial_istarget, A.eyelink_start_latency, A.eyelink_end_latency] = deal(nan(size(A.saccade_start_time)));
            for i=1:numel(trial_start_times)
                iThisTrial = find(A.saccade_end_time>trial_start_times(i) & A.saccade_start_time<trial_end_times(i));
                A.eyelink_trial_number(iThisTrial) = i;
                A.eyelink_trial_time(iThisTrial) = trial_start_times(i);
                A.eyelink_trial_istarget(iThisTrial) = object_istarget(i);
                A.eyelink_start_latency(iThisTrial) = A.saccade_start_time(iThisTrial)-trial_start_times(i);
                A.eyelink_end_latency(iThisTrial) = A.saccade_end_time(iThisTrial)-trial_start_times(i);
            end
            
            
        case 'sq' % squares dataset
            % EXTRACT
            A.saccade_start_time = x.saccade.start_time;
            A.saccade_end_time = x.saccade.end_time;
            A.saccade_start_x = x.saccade.start_position(:,1);
            A.saccade_start_y = x.saccade.start_position(:,2);
            A.saccade_end_x = x.saccade.end_position(:,1);
            A.saccade_end_y = x.saccade.end_position(:,2);
            % nSaccades = numel(A.saccade_start_time);

            % SORT
            % which saccades are in which epochs? what are their latencies in that epoch?
            [A.EEGset, A.EEGepoch, A.saccade_start_latency, A.saccade_end_latency] = deal(nan(size(A.saccade_start_time)));

            for j=1:numel(EEG)
                sacStartTimes = EyelinkToEegTimes(x.saccade.start_time,x)/x.eeg.eventsamplerate;
                [~, saccade_start_latency] = TimeToEpochNumber(EEG(j),sacStartTimes,sessionNumber,1);
                sacEndTimes = EyelinkToEegTimes(x.saccade.end_time,x)/x.eeg.eventsamplerate;
                [EEGepoch, saccade_end_latency] = TimeToEpochNumber(EEG(j),sacEndTimes,sessionNumber,1);

                A.saccade_start_latency(~isnan(EEGepoch)) = saccade_start_latency(~isnan(EEGepoch));
                A.saccade_end_latency(~isnan(EEGepoch)) = saccade_end_latency(~isnan(EEGepoch));
                A.EEGepoch(~isnan(EEGepoch)) = EEGepoch(~isnan(EEGepoch));
                A.EEGset(~isnan(EEGepoch)) = j;
                A.EEGsetnames{j} = EEG(j).setname;
            end

            % CALCULATE
            % duration
            A.saccade_duration = A.saccade_end_time - A.saccade_start_time;
            % distance
            x1 = A.saccade_start_x;
            x2 = A.saccade_end_x;
            y1 = A.saccade_start_y;
            y2 = A.saccade_end_y;
            A.saccade_length = sqrt((x2-x1).^2 + (y2-y1).^2);
            % average velocity
            A.saccade_velocity = A.saccade_length ./ A.saccade_duration;

            % LOCATE
            [~,~,~,A.saccade_end_disttoobj] = classify_squares_saccades(x);
            A.saccade_start_disttoobj = [Inf; A.saccade_end_disttoobj(1:end-1)]; % ONLY AN ESTIMATE!

            % TRACK    
            A.eyelink_session = repmat(x.session,size(A.saccade_start_time));
            trial_start_times = x.trial.start_time;
            trial_end_times = x.trial.end_time;
%             object_numbers = 1:length(trial_start_times);
            object_istarget = x.trial.is_target_trial;

            [A.eyelink_trial_number, A.eyelink_trial_time, A.eyelink_trial_istarget, A.eyelink_start_latency, A.eyelink_end_latency] = deal(nan(size(A.saccade_start_time)));
            for i=1:numel(trial_start_times)
                iThisTrial = find(A.saccade_end_time>trial_start_times(i) & A.saccade_start_time<trial_end_times(i));
                A.eyelink_trial_number(iThisTrial) = i;
                A.eyelink_trial_time(iThisTrial) = trial_start_times(i);
                A.eyelink_trial_istarget(iThisTrial) = object_istarget(i);
                A.eyelink_start_latency(iThisTrial) = A.saccade_start_time(iThisTrial)-trial_start_times(i);
                A.eyelink_end_latency(iThisTrial) = A.saccade_end_time(iThisTrial)-trial_start_times(i);
            end
    end
            
    
end

