function EEG = AddEeglabEvents(subject,session,input_suffix,events_rule,threshold,output_suffix)

% Removes the current events from the EEGLAB file and replaces them with
% new ones according to the events_rule specified in the code.  Then
% re-saves the data, overwriting the old file.
%
% AddEeglabEvents(subject,session,output_suffix,events_rule,threshold)
%
% - We think it's ok to overwrite the old file because you can always
% change the events back.  No data is lost during this script.
% - EEGLAB should already be started.  
%
% INPUTS:
% -subject and session are numbers that tell us which file to save
% -input_suffix is a string specifying the end of the filename to be loaded
% Full filename is '3DS-<subject>-<session><output_suffix>.set' [-filtered]
% -events_rule is a string indicating which set of events to use.  See
% MakeEventsMatrix for details. [EarlySaccades]
% -threshold is the other input to be passed to MakeEventsMatrix.  See
% that function's header for details. [200]
% -output_suffix is a string specifying the end of the filename to be saved
% Full filename is '3DS-<subject>-<session><output_suffix>.set' [input_suffix]
%
% OUTPUT:
% -EEG is the output file's eeglab struct.
%
% Created 8/2/10 by DJ.
% Updated 11/5/10 by DJ - switched to OddballTask events rule
% Updated 2/23/11 by DJ - made a function. 
% Updated 2/25/11 by DJ - added threshold. 
% Updated 3/1/11 by DJ - comments
% Updated 10/27/11 by DJ - added experimentType option to work with squares
%  experiments
% Updated 10/31/11 by DJ - separated input/output suffixes.
% Updated 2/28/13 by DJ - added SquaresFix compatibility.
% Updated 3/17/13 by DJ - added SquaresFix3 compatibility.

%% CHECK INPUTS AND SET UP
if nargin<3 || strcmp(input_suffix,'')
    input_suffix = '-filtered'; % follows '3DS-<subject>-<session>'
end
if nargin<4 || strcmp(events_rule,'')
    events_rule = 'EarlySaccades'; % see program MakeEventsMatrix for details
end
if nargin<5 || isempty(threshold)
    threshold = 200; % see program MakeEventsMatrix for details
end
if nargin<6 || strcmp(output_suffix,'')
    output_suffix = input_suffix;
end

if ~isempty(strfind(events_rule,'SquaresFix3')) % if the events rule contains the string 'Squares'
    experimentType = 'SquaresFix3';
elseif ~isempty(strfind(events_rule,'SquaresFix')) % if the events rule contains the string 'Squares'
    experimentType = 'SquaresFix';
elseif ~isempty(strfind(events_rule,'Squares')) % if the events rule contains the string 'Squares'
    experimentType = 'Squares';
else
    experimentType = '3DSearch';
end

% Set up
switch experimentType
    case '3DSearch'
        prefix = '3DS';        
    case 'Squares'
        prefix = 'sq';
    case 'SquaresFix'
        prefix = 'sf';
    case 'SquaresFix3'
        prefix = 'sf3';        
end
        
ALLEEG = [];

% Declare constants
data_dir = [cd '/'];
% data_dir = '/Users/dave/Documents/Data/3DSearch/';
filename = sprintf('%s-%d-%d%s.set',prefix,subject,session,input_suffix);
final_filename = sprintf('%s-%d-%d%s.set',prefix,subject,session,output_suffix);
% events_rule = 'OddballTask'; % see program MakeEventsMatrix for details
fprintf('---Adding events to %s according to rule %s\n',...
    filename,events_rule);
tic

%% Load data
EEG = pop_loadset('filename',filename,'filepath',data_dir);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );
load(sprintf('%s-%d-%d',prefix,subject,session)); % load trial structure x

%% Get events matrix and import into EEGLAB
events = MakeEventsMatrix(x,events_rule,threshold); % the times (in s) and codes of each event
assignin('base','events',events); % eeglab's importevent function grabs variable from base workspace
EEG = pop_importevent( EEG, 'append','no','event','events','fields',{'latency' 'type'},'timeunit',1,'optimalign','off');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

%% Save results
EEG = pop_saveset( EEG, 'filename',final_filename);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% Clean up
evalin('base','clear events'); % clean up
toc % Display elapsed time
% eeglab redraw
% clear data_dir events_rule x filename