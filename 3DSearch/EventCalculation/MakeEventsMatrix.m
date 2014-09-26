function events = MakeEventsMatrix(x,method,maxSeeTime)

% Create an events matrix that can be imported into EEGLAB and used to
% anchor epochs.
%
% events = MakeEventsMatrix(x,method,maxSeeTime)
% 
% INPUTS:
% - x is a 3DSearch data structure imported using ImportData.
% - method is a string saying which event types you would like.
%  Possible methods:
%   >'FiveSecondRule'
%       This method includes the first appearance of, and saccade to,
%       each object.  Objects that have been out of view for at least
%       five seconds are considered a new object.  Button press times and 
%       blink times are also included.
%   >'OddballTask'
%       This method includes ONLY THE FIRST appearance of and saccade to 
%       each object. Button press times, blink times, and ALL saccades 
%       are also included.
%   >'EarlySaccades' (requires maxSeeTime input)
%       This method is just like OddballTask, but it also includes SEES
%       events, which is the saccade end time (if the saccade was made
%       <maxSeeTime> ms after the stimulus) or the appear tie (otherwise).
%
% OUTPUTS:
% - events is an nx2 matrix, where n is the number of events.  The
% first column is the latency of each event, and the second column is the
% event code, which can be decoded by reading GetNumbers.m.
%
% Created 7/29/10 by DJ. 
% Updated 10/18/10 by DJ - added blink_times.
% Updated 11/5/10 by DJ - added OddballTask method
% Updated 11/22/10 by DJ - x.eeg.blink_times is now a column vector
% Updated 2/24/11 by DJ - added EarlySaccades and maxSeeTime input.
% Updated 6/2/11 by DJ - added SACCADE_START and SACCADE_END events
% Updated 7/28/11 by DJ - added BRAKE_BUTTON events
% Updated 10/27/11 by DJ - added Squares method
% Updated 11/16/11 by DJ - to match new squares data structures (x.sync)
% Updated 11/28/11 by DJ - added Squares-basic with text labels
% Updated 12/28/11 by DJ - added button.trialnum check to Squares-basic
% Updated 2/??/13 by DJ - added SquaresFix-basic rule
% Updated 5/17/13 by DJ - added SquaresFix3-basic rule

if nargin<3 % default
    maxSeeTime = 200;
end

switch method
    case 'FiveSecondRule'
        % This method includes the first appearance of, and saccade to,
        % each object.  Objects that have been out of view for at least
        % five seconds are considered a new object.  Button press times and 
        % blink times are also included.
        
        % Find the first saccade to each object
        % Include repeat objects that have been out of view for >5s
        [targetSac,distractorSac] = GetFirstSaccades(x,5);
        % Find the first time each object appeared in view
        % Include repeat objects that have been out of view for >5s
        [targetApp,distractorApp] = GetFirstAppears(x,5);
        times = [targetSac'; distractorSac'; targetApp'; distractorApp'; x.eeg.button_times'; x.eeg.blink_times'];
        % Convert to seconds
        times = times/x.eeg.eventsamplerate;
        % Find event codes
        GetNumbers;
        codes = [repmat(Numbers.SACCADE_TO+Numbers.TARGET,numel(targetSac),1); ...
            repmat(Numbers.SACCADE_TO+Numbers.DISTRACTOR,numel(distractorSac),1);...
            repmat(Numbers.ENTERS+Numbers.TARGET,numel(targetApp),1);...
            repmat(Numbers.ENTERS+Numbers.DISTRACTOR,numel(distractorApp),1);...
            repmat(Numbers.BUTTON,numel(x.eeg.button_times),1);
            repmat(Numbers.BLINK,numel(x.eeg.blink_times),1)];
        % sort times and events
        [times, order] = sort(times);
        codes = codes(order);
        % Construct events matrix
        events = [times, codes];
        
    case 'OddballTask'
        % This method includes ONLY THE FIRST appearance of and saccade to 
        % each object. Button press times, blink times, and ALL saccades 
        % are also included.
        
        % Find the first time each object appeared in view
        % Do not include repeat objects 
        [targetApp,distractorApp] = GetFirstAppears(x,Inf);
        % Find the first saccade to each object
        % Do not include repeat objects 
        [targetSac,distractorSac] = GetFirstSaccades(x,Inf);
        saccadeTimesEeg = EyelinkToEegTimes(x.eyelink.saccade_times,x);
        sacStartTimesEeg = EyelinkToEegTimes(x.eyelink.saccade_start_times,x);
        times = [targetSac'; distractorSac'; targetApp'; distractorApp'; ...
            x.eeg.button_times; x.eeg.brake_button_times; x.eeg.blink_times; saccadeTimesEeg; sacStartTimesEeg];
        % Convert to seconds
        times = times/x.eeg.eventsamplerate;
        % Find event codes
        GetNumbers;
        codes = [repmat(Numbers.SACCADE_TO+Numbers.TARGET,numel(targetSac),1); ...
            repmat(Numbers.SACCADE_TO+Numbers.DISTRACTOR,numel(distractorSac),1);...
            repmat(Numbers.ENTERS+Numbers.TARGET,numel(targetApp),1);...
            repmat(Numbers.ENTERS+Numbers.DISTRACTOR,numel(distractorApp),1);...
            repmat(Numbers.BUTTON,numel(x.eeg.button_times),1);...
            repmat(Numbers.BRAKE_BUTTON,numel(x.eeg.brake_button_times),1);...
            repmat(Numbers.BLINK,numel(x.eeg.blink_times),1);...
            repmat(Numbers.SACCADE_END,numel(saccadeTimesEeg),1);...
            repmat(Numbers.SACCADE_START,numel(sacStartTimesEeg),1)];
            
        % sort times and events
        [times, order] = sort(times);
        codes = codes(order);
        % Construct events matrix
        events = [times, codes];
    
    case 'EarlySaccades'
        % This method includes ONLYT THE FIRST appearance of and saccade to 
        % each object. Button press times, blink times, ALL saccades, 
        % simulated button presses (for distractor trials), and "see" 
        % times (saccade if within <maxSeeTime> ms of target appearance, 
        % stim time otherwise) are also included.
        
        % Find the first time each object appeared in view
        % Do not include repeat objects 
        [targetApp,distractorApp] = GetFirstAppears(x,Inf);
        % Find the first saccade to each object
        % Do not include repeat objects 
        [targetSac,distractorSac] = GetFirstSaccades(x,Inf);
        % If a saccade was made <maxSeeTime> ms after appear, pick saccade time.  
        % Otherwise, pick appear time.
        targetTime = zeros(1,numel(targetApp));
        for i=1:numel(targetApp)
            iInRange = find(targetSac > targetApp(i) & (targetSac-targetApp(i)) < (maxSeeTime*x.eeg.eventsamplerate/1000)); % adjust to EEG samples
            if isempty(iInRange)
                targetTime(i) = targetApp(i);
            else
                targetTime(i) = targetSac(iInRange);
            end
        end
        distractorTime = zeros(1,numel(distractorApp));
        for i=1:numel(distractorApp)
            iInRange = find(distractorSac > distractorApp(i) & (distractorSac-distractorApp(i)) < (maxSeeTime*x.eeg.eventsamplerate/1000)); % adjust to EEG samples
            if isempty(iInRange)
                distractorTime(i) = distractorApp(i);
            else
                distractorTime(i) = distractorSac(iInRange);
            end
        end
        % Get distribution of reaction times, and use to simulate button
        % presses for distractor trials.
        RT = [];
        for i=1:numel(targetApp)
            iInRange = find(x.eeg.button_times > targetApp(i) & (x.eeg.button_times-targetApp(i)) < (2*x.eeg.eventsamplerate)); % 2 seconds, adjusted to EEG samples
            if ~isempty(iInRange)                
                RT = [RT x.eeg.button_times(iInRange(1))-targetApp(i)];                
            end
        end
        % Simulate button presses by selecting a random reaction time (with
        % replacement).
        if isempty(RT) % if there were no button presses
           simButtonTimes = [];
        else % if there were button presses
            simButtonTimes = zeros(1,numel(distractorApp));
            for i=1:numel(distractorApp)              
                simButtonTimes(i) = distractorApp(i)+RT(randi(numel(RT)));
            end      
        end
        
        saccadeTimesEeg = EyelinkToEegTimes(x.eyelink.saccade_times,x);
        sacStartTimesEeg = EyelinkToEegTimes(x.eyelink.saccade_start_times,x);
        times = [targetSac'; distractorSac'; targetApp'; distractorApp'; ...
            targetTime'; distractorTime'; x.eeg.button_times; simButtonTimes'; ...
            x.eeg.brake_button_times; x.eeg.blink_times; saccadeTimesEeg; sacStartTimesEeg];
        % Convert to seconds
        times = times/x.eeg.eventsamplerate;
        % Find event codes
        GetNumbers;
        codes = [repmat(Numbers.SACCADE_TO+Numbers.TARGET,numel(targetSac),1); ...
            repmat(Numbers.SACCADE_TO+Numbers.DISTRACTOR,numel(distractorSac),1);...
            repmat(Numbers.ENTERS+Numbers.TARGET,numel(targetApp),1);...
            repmat(Numbers.ENTERS+Numbers.DISTRACTOR,numel(distractorApp),1);...
            repmat(Numbers.SEES+Numbers.TARGET,numel(targetTime),1);...
            repmat(Numbers.SEES+Numbers.DISTRACTOR,numel(distractorTime),1);...            
            repmat(Numbers.BUTTON,numel(x.eeg.button_times),1);...
            repmat(Numbers.SIM_BUTTON,numel(simButtonTimes),1);...
            repmat(Numbers.BRAKE_BUTTON,numel(x.eeg.brake_button_times),1);...
            repmat(Numbers.BLINK,numel(x.eeg.blink_times),1);...
            repmat(Numbers.SACCADE_END,numel(saccadeTimesEeg),1);...
            repmat(Numbers.SACCADE_START,numel(sacStartTimesEeg),1)];
        % sort times and events
        [times, order] = sort(times);
        codes = codes(order);
        % Construct events matrix
        events = [times, codes];
        
    case 'Squares'
        % This method includes appearances, classified saccades, button 
        % presses, and blinks.        
        
        Constants = GetSquaresConstants;        
        % Get saccade codes
        sacCodes = x.saccade.class;
        sacCodes(isnan(sacCodes)) = Constants.OTHER;
        % Get fixation codes
        fixDir = repmat(Constants.OTHERSIDE,size(x.fixation.start_time));
        fixDir(sqrt((x.fixation.position(:,1)-Constants.LEFTCROSS_X).^2+(x.fixation.position(:,2)-Constants.LEFTCROSS_Y).^2) < x.pixel_threshold) = Constants.LEFTSIDE;
        fixDir(sqrt((x.fixation.position(:,1)-Constants.RIGHTCROSS_X).^2+(x.fixation.position(:,2)-Constants.RIGHTCROSS_Y).^2) < x.pixel_threshold) = Constants.RIGHTSIDE;
        % Construct events matrix
        elEvents = [x.block_start_time,Constants.BLOCK_START;...
            x.trial.fix_time,Constants.FIXCROSS_ON+x.trial.is_target_trial;...
            x.trial.start_time,Constants.TRIAL_START+x.trial.is_target_trial;...
            x.trial.end_time,Constants.TRIAL_END+x.trial.is_target_trial;...
            x.saccade.start_time,Constants.SACCADESTART_BASE+sacCodes;...
            x.saccade.end_time,Constants.SACCADEEND_BASE+sacCodes;...
            x.button.time,Constants.BUTTON_BASE+x.button.number;...
            x.fixation.start_time,Constants.FIXSTART_BASE+fixDir;...
            x.fixation.end_time,Constants.FIXEND_BASE+fixDir;...
            x.blink.start_time,repmat(Constants.BLINKSTART,size(x.blink.start_time));...
            x.blink.end_time,repmat(Constants.BLINKEND,size(x.blink.end_time))];
          
        % separate times and codes
        eegTimes = EyelinkToEegTimes(elEvents(:,1),x.sync.eyelink,x.sync.eeg);
        codes = elEvents(:,2);
        % Convert to seconds
        times = eegTimes/x.eeg.eventsamplerate;           
        
        % sort times and events
        [times, order] = sort(times);
        codes = codes(order);
        % Construct events matrix
        events = [times, codes];
        
    case 'Squares-eog'
        Constants = GetSquaresConstants;
        fixDir = repmat(Constants.OTHERSIDE,size(x.fixation.start_time));
        fixDir(sqrt((x.fixation.position(:,1)-Constants.LEFTCROSS_X).^2+(x.fixation.position(:,2)-Constants.LEFTCROSS_Y).^2) < x.pixel_threshold) = Constants.LEFTSIDE;
        fixDir(sqrt((x.fixation.position(:,1)-Constants.RIGHTCROSS_X).^2+(x.fixation.position(:,2)-Constants.RIGHTCROSS_Y).^2) < x.pixel_threshold) = Constants.RIGHTSIDE;
        
        elEvents = [x.fixation.start_time,Constants.FIXSTART_BASE+fixDir;...
            x.fixation.end_time,Constants.FIXEND_BASE+fixDir;...
            x.blink.start_time,repmat(Constants.BLINKSTART,size(x.blink.start_time));...
            x.blink.end_time,repmat(Constants.BLINKEND,size(x.blink.end_time))];
        % separate times and codes
        eegTimes = EyelinkToEegTimes(elEvents(:,1),x.sync.eyelink,x.sync.eeg);
        codes = elEvents(:,2);
        % Convert to seconds
        times = eegTimes/x.eeg.eventsamplerate;           
        
        % sort times and events
        [times, order] = sort(times);
        codes = codes(order);
        % Construct events matrix
        events = [times, codes];
        
    case 'Squares-basic'
        % This method includes appearances, basic saccades, button 
        % presses, and blinks.        
        
        Constants = GetSquaresConstants;        
        % Get fixation codes
        fixDir = nan(size(x.fixation.start_time));
        fixDir(sqrt((x.fixation.position(:,1)-Constants.LEFTCROSS_X).^2+(x.fixation.position(:,2)-Constants.LEFTCROSS_Y).^2) < x.pixel_threshold) = Constants.LEFTSIDE;
        fixDir(sqrt((x.fixation.position(:,1)-Constants.RIGHTCROSS_X).^2+(x.fixation.position(:,2)-Constants.RIGHTCROSS_Y).^2) < x.pixel_threshold) = Constants.RIGHTSIDE;
        fix_start = x.fixation.start_time(~isnan(fixDir));
        fix_end = x.fixation.end_time(~isnan(fixDir));
        fixDir = fixDir(~isnan(fixDir));
        % Get saccade codes - ALL ZEROS!
        sacCodes = zeros(size(x.saccade.class));  
        % Get button codes
        if ~isfield(x.button,'trialnum')
            warning('x.button.trialnum not found... some non-response button presses may be incorrectly labeled!');
            x.button.trialnum = zeros(size(x.button.number));
        end
        buttonCodes = repmat(2,size(x.button.time)); % start with everything labeled as "other button"
        buttonCodes(~isnan(x.button.trialnum) & x.button.number==Constants.NONTARGET_BUTTON) = 0;
        buttonCodes(~isnan(x.button.trialnum) & x.button.number==Constants.TARGET_BUTTON) = 1;
        % Construct events matrix
        elEvents = [x.block_start_time,Constants.BLOCK_START;...
            x.trial.fix_time,Constants.FIXCROSS_ON+x.trial.is_target_trial;...
            x.trial.start_time,Constants.TRIAL_START+x.trial.is_target_trial;...
            x.trial.end_time,Constants.TRIAL_END+x.trial.is_target_trial;...
            x.saccade.start_time,Constants.SACCADESTART_BASE+sacCodes;...
            x.saccade.end_time,Constants.SACCADEEND_BASE+sacCodes;...
            x.button.time,Constants.BUTTON_BASE+buttonCodes;...
            x.blink.start_time,repmat(Constants.BLINKSTART,size(x.blink.start_time));...
            x.blink.end_time,repmat(Constants.BLINKEND,size(x.blink.end_time));...
            fix_start,Constants.FIXSTART_BASE+fixDir;...
            fix_end,Constants.FIXEND_BASE+fixDir];
          
        % separate times and codes
        eegTimes = EyelinkToEegTimes(elEvents(:,1),x.sync.eyelink,x.sync.eeg);
        codes = Constants.EVENTNAMES(elEvents(:,2))'; % text codes
        if any(cellfun('isempty',codes))
            warning('DaveJ:eventnotfound','Some event codes were not recognized!')
        end
        % Convert to seconds
        times = eegTimes/x.eeg.eventsamplerate;           
        
        % sort times and events
        [times, order] = sort(times);
        codes = codes(order);
        % Construct events matrix
        events = [num2cell(times), codes];
    
    
    case {'SquaresFix-basic', 'SquaresFix3-basic'}
        % This method includes appearances, button presses, and blinks.        
        
        Constants = GetSquaresConstants;        
        % Get saccade codes - ALL ZEROS!
        sacCodes = zeros(size(x.saccade.start_time));  
        % Get button codes
        if ~isfield(x.button,'trialnum')
            warning('x.button.trialnum not found... some non-response button presses may be incorrectly labeled!');
            x.button.trialnum = zeros(size(x.button.number));
        end
        buttonCodes = repmat(2,size(x.button.time)); % start with everything labeled as "other button"
        buttonCodes(~isnan(x.button.trialnum) & x.button.number==Constants.NONTARGET_BUTTON) = 0;
        buttonCodes(~isnan(x.button.trialnum) & x.button.number==Constants.TARGET_BUTTON) = 1;
        % Construct events matrix
        elEvents = [x.block_start_time,Constants.BLOCK_START;...
            x.trial.fix_time,Constants.FIXCROSS_ON+x.trial.is_target_trial;...
            x.trial.start_time,Constants.TRIAL_START+x.trial.is_target_trial;...
            x.trial.end_time,Constants.TRIAL_END+x.trial.is_target_trial;...
            x.saccade.start_time,Constants.SACCADESTART_BASE+sacCodes;...
            x.saccade.end_time,Constants.SACCADEEND_BASE+sacCodes;...
            x.button.time,Constants.BUTTON_BASE+buttonCodes;...
            x.blink.start_time,repmat(Constants.BLINKSTART,size(x.blink.start_time));...
            x.blink.end_time,repmat(Constants.BLINKEND,size(x.blink.end_time));
            x.trial.square_time(:),repmat(Constants.SQUARE_ON,size(x.trial.square_time(:)))];
          
        % separate times and codes
        eegTimes = EyelinkToEegTimes(elEvents(:,1),x.sync.eyelink,x.sync.eeg);
        codes = Constants.EVENTNAMES(elEvents(:,2))'; % text codes
        if any(cellfun('isempty',codes))
            warning('DaveJ:eventnotfound','Some event codes were not recognized!')
        end
        % Convert to seconds
        times = eegTimes/x.eeg.eventsamplerate;           
        
        % sort times and events
        [times, order] = sort(times);
        codes = codes(order);
        % Construct events matrix
        events = [num2cell(times), codes];

    otherwise
        error('Method not recognized!');
end