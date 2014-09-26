function varargout = EpochSquaresGui(varargin)
% EPOCHSQUARESGUI M-file for EpochSquaresGui.fig
%      EPOCHSQUARESGUI, by itself, creates a new EPOCHSQUARESGUI or raises the existing
%      singleton*.
%
%      H = EPOCHSQUARESGUI returns the handle to a new EPOCHSQUARESGUI or the handle to
%      the existing singleton*.
%
%      EPOCHSQUARESGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EPOCHSQUARESGUI.M with the given input arguments.
%
%      EPOCHSQUARESGUI('Property','Value',...) creates a new EPOCHSQUARESGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EpochSquaresGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EpochSquaresGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EpochSquaresGui

% Last Modified by GUIDE v2.5 28-Nov-2011 19:50:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EpochSquaresGui_OpeningFcn, ...
                   'gui_OutputFcn',  @EpochSquaresGui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before EpochSquaresGui is made visible.
function EpochSquaresGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EpochSquaresGui (see VARARGIN)

% Choose default command line output for EpochSquaresGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EpochSquaresGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EpochSquaresGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        BUTTON FUNCTIONS                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in button_import.
function button_import_Callback(hObject, eventdata, handles)
% hObject    handle to button_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get info from GUI fields
subject_str = get(handles.edit_subject,'string');
handles.subject = str2num(subject_str);
sessions_str = get(handles.edit_sessions,'string');
handles.sessions = str2num(sessions_str);
eegsessions_str = get(handles.edit_eegsessions,'string');
handles.eegsessions = str2num(eegsessions_str);
filetype_str = get(handles.popup_filetype,'string');
handles.filetype = filetype_str{get(handles.popup_filetype,'value')};
eventsrule_str = get(handles.popup_eventsrule,'string');
handles.eventsrule = eventsrule_str{get(handles.popup_eventsrule,'value')};
pixelthresh_str = get(handles.edit_pixelthresh,'string');
handles.pixelthresh = str2num(pixelthresh_str);
handles.singlesuffix = get(handles.edit_singlesuffix,'string');
handles.combosuffix = get(handles.edit_combosuffix,'string');
guidata(hObject,handles); % Save variables


for i=1:numel(handles.sessions)
    % Run functions as in ImportData
    Import_squares_data(handles.subject,handles.sessions(i),handles.eegsessions(i),handles.filetype,handles.pixelthresh);
    ImportToEeglab(handles.subject,handles.sessions(i),handles.eegsessions(i),handles.filetype,'Squares');
    AddEeglabEvents(handles.subject,handles.sessions(i),handles.singlesuffix,handles.eventsrule);
end
%Combine sessions
CombineEeglabSessions(handles.subject,handles.sessions,handles.singlesuffix,handles.combosuffix,'Squares');



% --- Executes on button press in button_recalcevents.
function button_recalcevents_Callback(hObject, eventdata, handles)
% hObject    handle to button_recalcevents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get info from GUI fields
subject_str = get(handles.edit_subject,'string');
handles.subject = str2num(subject_str);
sessions_str = get(handles.edit_sessions,'string');
handles.sessions = str2num(sessions_str);
eventsrule_str = get(handles.popup_eventsrule,'string');
handles.eventsrule = eventsrule_str{get(handles.popup_eventsrule,'value')};
pixelthresh_str = get(handles.edit_pixelthresh,'string');
handles.pixelthresh = str2num(pixelthresh_str);
handles.singlesuffix = get(handles.edit_singlesuffix,'string');
handles.combosuffix = get(handles.edit_combosuffix,'string');
guidata(hObject,handles); % Save variables

% Use info to change files!
ChangeSacToObjThresholds_squares(handles.subject,handles.sessions,handles.pixelthresh,...
    handles.eventsrule,handles.singlesuffix,handles.combosuffix);



% --- Executes on button press in button_removeduds.
function button_removeduds_Callback(hObject, eventdata, handles)
% hObject    handle to button_removeduds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get info from GUI fields
subject_str = get(handles.edit_subject,'string');
handles.subject = str2num(subject_str);
duds_str = get(handles.edit_duds,'string');
eval(sprintf('handles.duds = {%s};',duds_str)); % make duds a cell array of strings
handles.inputsuffix = get(handles.edit_inputsuffix,'string');
handles.outputsuffix = get(handles.edit_outputsuffix,'string');
guidata(hObject,handles); % Save variables

% Remove bad electrodes
RemoveElectrodes(sprintf('sq-%d%s.set',handles.subject,handles.inputsuffix),sprintf('sq-%d%s.set',handles.subject,handles.outputsuffix),handles.duds);



% --- Executes on button press in button_removefrontal.
function button_removefrontal_Callback(hObject, eventdata, handles)
% hObject    handle to button_removefrontal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get info from GUI fields
subject_str = get(handles.edit_subject,'string');
handles.subject = str2num(subject_str);
handles.inputsuffix = get(handles.edit_inputsuffix,'string');
handles.outputsuffix = get(handles.edit_outputsuffix,'string');
filetype_str = get(handles.popup_filetype,'string');
handles.filetype = filetype_str{get(handles.popup_filetype,'value')};
guidata(hObject,handles); % Save variables

% Remove bad electrodes
RemoveFrontalElectrodes(sprintf('sq-%d%s.set',handles.subject,handles.inputsuffix),sprintf('sq-%d%s.set',handles.subject,handles.outputsuffix),handles.filetype);


% --- Executes on button press in button_removeblinks.
function button_removeblinks_Callback(hObject, eventdata, handles)
% hObject    handle to button_removeblinks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get info from GUI fields
subject_str = get(handles.edit_subject,'string');
handles.subject = str2num(subject_str);
handles.inputsuffix = get(handles.edit_inputsuffix,'string');
handles.outputsuffix = get(handles.edit_outputsuffix,'string');
filetype_str = get(handles.popup_filetype,'string');
handles.filetype = filetype_str{get(handles.popup_filetype,'value')};
guidata(hObject,handles); % Save variables

% Remove blink artifacts
blink_expand = [-50 100]; % ms to expand blinks outside the eyelink saccade events
outside_window = Inf;%25; % ms to average outside the blink window for mean subtraction
RemoveBlinkArtifacts(sprintf('sq-%d%s.set',handles.subject,handles.inputsuffix),...
    sprintf('sq-%d%s.set',handles.subject,handles.outputsuffix),blink_expand,outside_window);


% --- Executes on button press in button_removeeog.
function button_removeeog_Callback(hObject, eventdata, handles)
% hObject    handle to button_removeeog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get info from GUI fields
subject_str = get(handles.edit_subject,'string');
handles.subject = str2num(subject_str);
handles.inputsuffix = get(handles.edit_inputsuffix,'string');
handles.outputsuffix = get(handles.edit_outputsuffix,'string');
filetype_str = get(handles.popup_filetype,'string');
handles.filetype = filetype_str{get(handles.popup_filetype,'value')};
guidata(hObject,handles); % Save variables

% Remove blink artifacts
fix_expand = [0 0]; % ms to expand blinks outside the eyelink saccade events
outside_window = Inf; % ms to average outside the blink window for mean subtraction
RemoveEogArtifacts(sprintf('sq-%d%s.set',handles.subject,handles.inputsuffix),...
    sprintf('sq-%d%s.set',handles.subject,handles.outputsuffix),fix_expand,outside_window);


% --- Executes on button press in button_epoch.
function button_epoch_Callback(hObject, eventdata, handles)
% hObject    handle to button_epoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get info from GUI fields
handles.subject = str2double(get(handles.edit_subject,'string'));
handles.version = get(handles.edit_version,'string');
handles.eegFilename = sprintf('sq-%d-%s.set',handles.subject,handles.version);
reference_str = get(handles.popup_reference,'string');
handles.reference = reference_str{get(handles.popup_reference,'value')};
if strcmp(handles.reference,'Select Reference') % use default
    handles.reference = 'None';
end
epochtimes_str = get(handles.edit_epochtimes,'string');
handles.epochtimes = str2num(epochtimes_str);
baselinetimes_str = get(handles.edit_baselinetimes,'string');
handles.baselinetimes = str2num(baselinetimes_str);
baseevent_str = get(handles.popup_baseevent,'string');
handles.baseevent = baseevent_str{get(handles.popup_baseevent,'value')};
handles.removeIncorrect = get(handles.check_incorrect,'value');
handles.removeBlinks = get(handles.check_blinks,'value');
handles.useCutoff = get(handles.check_cutoff,'value');
cutoff_str = get(handles.edit_cutoff,'string');
handles.cutoff = str2num(cutoff_str);
tstimon_str = get(handles.edit_t_win_stimon,'string');
handles.tstimon = str2num(tstimon_str);
tepoch_str = get(handles.edit_t_win_epoch,'string');
handles.tepoch = str2num(tepoch_str);

% Load EEG data
fprintf('Loading behavior for S%d...\n',handles.subject);
y = loadBehaviorData(handles.subject); 
fprintf('Loading file %s...\n',handles.eegFilename);
EEG = pop_loadset(handles.eegFilename);

% Get events from GUI checkboxes
handles.eventnames = {};
Constants = GetSquaresConstants; % define numbers struct containing Squares constants
checknames = {'trialstartd','trialstartt','dist0t','dist1t','integ','compl','buttond','buttont'};
eventnames = {'TrialStart-D','TrialStart-T','Dist-0T','Dist-1T','Integ','Compl','Button-D','Button-T'};
for i=1:numel(checknames)
    if get(handles.(['check_event_' checknames{i}]),'value')
        handles.eventnames = [handles.eventnames, eventnames(i)];        
        [elTimes,codes] = UseEventRule(y,eventnames{i});
        EEG = AddEeglabEvents_MultiSession(EEG,y,elTimes,codes);
    end
end
guidata(hObject,handles); % Save variables



% Run Epoching - because we're using EEGLAB, we must evaluate these
% expressions in the base workspace, where EEGLAB is running!

% USE BASELINE RELATIVE TO STIMULUS APPEARANCE
if strcmp(handles.baseevent,'Stimulus')
    % epoch to appear time WITH baseline subtraction and wide epoch times
    [ALLEEG] = EpochEeglabFile(EEG,[-2000 5000],handles.baselinetimes,...
        {'TrialStart-D', 'TrialStart-T'},...
        {'baseline=appear'},handles.reference); 
    %         [Constants.TRIAL_START,
    %         Constants.TRIAL_START+Constants.ISTARGET],...
    EEG = ALLEEG(end);
    % Remove incorrect trials, if requested
    if handles.removeIncorrect
        EEG = RemoveIncorrectTrials(EEG,'Squares');
    end
    if handles.removeBlinks && diff(handles.tstimon)>0
        EEG = RemoveEyeBlinkTrials(EEG,Constants.BLINKSTART,handles.tstimon,[0 -1],[0 -1],'Squares');        
    end
    if handles.useCutoff && diff(handles.tstimon)>0
        EEG = ExcludeTrialsByCutoff(EEG,handles.cutoff,handles.tstimon,[0 -1],[0 -1],'Squares');
    end
    % re-epoch to individual events with desired epoch times and WITHOUT baseline subtraction
    [ALLEEG, EEG, CURRENTSET] = EpochEeglabFile(EEG,handles.epochtimes,[],...
        handles.eventnames,handles.eventnames,handles.reference); 

% OR USE BASELINE RELATIVE TO EACH INDIVIDUAL EVENT
else
    % epoch to appear time WITHOUT baseline subtraction and wide epoch times
    [ALLEEG] = EpochEeglabFile(EEG,[-2000 5000],[],...
        {'TrialStart-D', 'TrialStart-T'},...
        {'baseline=appear'},handles.reference);
    %         [Constants.TRIAL_START,
    %         Constants.TRIAL_START+Constants.ISTARGET],...
    EEG = ALLEEG(end);
    % Remove incorrect trials, if requested
    if handles.removeIncorrect
        EEG = RemoveIncorrectTrials(EEG,'Squares');
    end
    if handles.removeBlinks && diff(handles.tstimon)>0
        EEG = RemoveEyeBlinkTrials(EEG,Constants.BLINKSTART,handles.tstimon,[0 -1],[0 -1],'Squares');        
    end
    if handles.useCutoff && diff(handles.tstimon)>0
        EEG = ExcludeTrialsByCutoff(EEG,handles.cutoff,handles.tstimon,[0 -1],[0 -1],'Squares');
    end
    % re-epoch to individual events with desired epoch times and WITH baseline subtraction
    [ALLEEG, EEG, CURRENTSET] = EpochEeglabFile(EEG,handles.epochtimes,handles.baselinetimes,handles.eventnames,handles.eventnames,handles.reference);
end
    
% PERFORM EPOCHING
nFiles = length(ALLEEG);
for i=2:nFiles
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',i,'study',0);
    fprintf('------Cleaning up %s dataset...\n',handles.eventnames{i-1});
    if handles.removeBlinks && diff(handles.tepoch)>0
        EEG = RemoveEyeBlinkTrials(EEG,Constants.BLINKSTART,handles.tepoch,[0 -1],[0 -1],'Squares');        
    end
    if handles.useCutoff && diff(handles.tepoch)>0
        EEG = ExcludeTrialsByCutoff(EEG,handles.cutoff,handles.tepoch,[0 -1],[0 -1],'Squares');
    end
    %  Store results (but do not save)
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off'); 
    EEG = eeg_checkset( EEG );
end
% Send results to base workspace so user can explore (and use them for
% analysis parts of this gui)
assignin('base','ALLEEG',ALLEEG);
assignin('base','EEG',EEG);
assignin('base','CURRENTSET',CURRENTSET);
evalin('base','eeglab redraw');

% Send names to popups
set(handles.popup_set1,'string',[{'Select Dataset 1'},handles.eventnames],'value',1);
set(handles.popup_set2,'string',[{'Select Dataset 2'},handles.eventnames],'value',1);
disp('------Done!');

% --- Executes on button press in button_ploterps.
function button_ploterps_Callback(hObject, eventdata, handles)
% hObject    handle to button_ploterps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get info from GUI fields
handles.set1 = get(handles.popup_set1,'value');
handles.set2 = get(handles.popup_set2,'value');
guidata(hObject,handles); % Save variables

% Draw ERP
ALLEEG = evalin('base','ALLEEG'); % Get ALLEEG from base workspace
MakeErpMovie(ALLEEG,handles.set1,handles.set2); % Make ERPs


% --- Executes on button press in button_run.
function button_run_Callback(hObject, eventdata, handles)
% hObject    handle to button_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get info from GUI fields
handles.set1 = get(handles.popup_set1,'value');
handles.set2 = get(handles.popup_set2,'value');
handles.doloo = get(handles.check_doloo,'value');
handles.bootstrap = get(handles.check_bootstrap,'value');
spw_str = get(handles.edit_spw,'string'); % Samples per window
handles.spw = str2num(spw_str);
sps_str = get(handles.edit_sps,'string'); % Samples per shift
handles.sps = str2num(sps_str);
guidata(hObject,handles); % Save variables

% Select which datasets to use (distractors, then targets for proper forward model).  
GetNumbers;

% Run Analysis - because we're using EEGLAB, we must evaluate these
% expressions in the base workspace, where EEGLAB is running!
ALLEEG = evalin('base','ALLEEG'); % Get ALLEEG from base workspace
[~, time, Azloo] = MakeLogRegMovie(ALLEEG,[handles.set1, handles.set2],1:ALLEEG(handles.set1).nbchan,handles.spw,handles.sps,handles.doloo,handles.bootstrap);

if handles.doloo
    % Save the results in the 'LOO.mat' file
    SaveLooResults(ALLEEG(handles.set1),ALLEEG(handles.set2),time,Azloo,handles.spw,handles.sps,datestr(now),handles.bootstrap);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         EDIT FUNCTIONS                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_subject_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subject as text
%        str2double(get(hObject,'String')) returns contents of edit_subject as a double

% Get subject number
subject_str = get(handles.edit_subject,'string');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        UNUSED FUNCTIONS                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function edit_subject_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sessions_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sessions as text
%        str2double(get(hObject,'String')) returns contents of edit_sessions as a double


% --- Executes during object creation, after setting all properties.
function edit_sessions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_filetype.
function popup_filetype_Callback(hObject, eventdata, handles)
% hObject    handle to popup_filetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_filetype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_filetype


% --- Executes during object creation, after setting all properties.
function popup_filetype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_filetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function edit_eyesession_Callback(hObject, eventdata, handles)
% hObject    handle to edit_eyesession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_eyesession as text
%        str2double(get(hObject,'String')) returns contents of edit_eyesession as a double


% --- Executes during object creation, after setting all properties.
function edit_eyesession_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_eyesession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_blinks.
function check_blinks_Callback(hObject, eventdata, handles)
% hObject    handle to check_blinks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_blinks


% --- Executes on button press in check_cutoff.
function check_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to check_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_cutoff



function edit_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cutoff as text
%        str2double(get(hObject,'String')) returns contents of edit_cutoff as a double


% --- Executes during object creation, after setting all properties.
function edit_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_set1.
function popup_set1_Callback(hObject, eventdata, handles)
% hObject    handle to popup_set1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_set1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_set1


% --- Executes during object creation, after setting all properties.
function popup_set1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_set1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_doloo.
function check_doloo_Callback(hObject, eventdata, handles)
% hObject    handle to check_doloo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_doloo


% --- Executes on button press in check_bootstrap.
function check_bootstrap_Callback(hObject, eventdata, handles)
% hObject    handle to check_bootstrap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_bootstrap



function edit_spw_Callback(hObject, eventdata, handles)
% hObject    handle to edit_spw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_spw as text
%        str2double(get(hObject,'String')) returns contents of edit_spw as a double


% --- Executes during object creation, after setting all properties.
function edit_spw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_spw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sps_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sps as text
%        str2double(get(hObject,'String')) returns contents of edit_sps as a double


% --- Executes during object creation, after setting all properties.
function edit_sps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_eegsessions_Callback(hObject, eventdata, handles)
% hObject    handle to edit_eegsessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_eegsessions as text
%        str2double(get(hObject,'String')) returns contents of edit_eegsessions as a double


% --- Executes during object creation, after setting all properties.
function edit_eegsessions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_eegsessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_version.
function popup_version_Callback(hObject, eventdata, handles)
% hObject    handle to popup_version (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_version contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_version


% --- Executes during object creation, after setting all properties.
function popup_version_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_version (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_reference.
function popup_reference_Callback(hObject, eventdata, handles)
% hObject    handle to popup_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_reference contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_reference


% --- Executes during object creation, after setting all properties.
function popup_reference_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_eventsrule.
function popup_eventsrule_Callback(hObject, eventdata, handles)
% hObject    handle to popup_eventsrule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_eventsrule contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_eventsrule


% --- Executes during object creation, after setting all properties.
function popup_eventsrule_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_eventsrule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_singlesuffix_Callback(hObject, eventdata, handles)
% hObject    handle to edit_singlesuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_singlesuffix as text
%        str2double(get(hObject,'String')) returns contents of edit_singlesuffix as a double


% --- Executes during object creation, after setting all properties.
function edit_singlesuffix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_singlesuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_combosuffix_Callback(hObject, eventdata, handles)
% hObject    handle to edit_combosuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_combosuffix as text
%        str2double(get(hObject,'String')) returns contents of edit_combosuffix as a double


% --- Executes during object creation, after setting all properties.
function edit_combosuffix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_combosuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pixelthresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pixelthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pixelthresh as text
%        str2double(get(hObject,'String')) returns contents of edit_pixelthresh as a double


% --- Executes during object creation, after setting all properties.
function edit_pixelthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pixelthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_timelimits_Callback(hObject, eventdata, handles)
% hObject    handle to edit_timelimits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_timelimits as text
%        str2double(get(hObject,'String')) returns contents of edit_timelimits as a double


% --- Executes during object creation, after setting all properties.
function edit_timelimits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_timelimits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_duds_Callback(hObject, eventdata, handles)
% hObject    handle to edit_duds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_duds as text
%        str2double(get(hObject,'String')) returns contents of edit_duds as a double


% --- Executes during object creation, after setting all properties.
function edit_duds_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_duds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_inputsuffix_Callback(hObject, eventdata, handles)
% hObject    handle to edit_inputsuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_inputsuffix as text
%        str2double(get(hObject,'String')) returns contents of edit_inputsuffix as a double


% --- Executes during object creation, after setting all properties.
function edit_inputsuffix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_inputsuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_outputsuffix_Callback(hObject, eventdata, handles)
% hObject    handle to edit_outputsuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_outputsuffix as text
%        str2double(get(hObject,'String')) returns contents of edit_outputsuffix as a double


% --- Executes during object creation, after setting all properties.
function edit_outputsuffix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_outputsuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxseetime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxseetime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxseetime as text
%        str2double(get(hObject,'String')) returns contents of edit_maxseetime as a double


% --- Executes during object creation, after setting all properties.
function edit_maxseetime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxseetime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_epochtimes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_epochtimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_epochtimes as text
%        str2double(get(hObject,'String')) returns contents of edit_epochtimes as a double


% --- Executes during object creation, after setting all properties.
function edit_epochtimes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_epochtimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_baselinetimes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_baselinetimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_baselinetimes as text
%        str2double(get(hObject,'String')) returns contents of edit_baselinetimes as a double


% --- Executes during object creation, after setting all properties.
function edit_baselinetimes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_baselinetimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_event_integ.
function check_event_integ_Callback(hObject, eventdata, handles)
% hObject    handle to check_event_integ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_event_integ


% --- Executes on button press in check_event_dist0t.
function check_event_dist0t_Callback(hObject, eventdata, handles)
% hObject    handle to check_event_dist0t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_event_dist0t


% --- Executes on button press in check_event_compl.
function check_event_compl_Callback(hObject, eventdata, handles)
% hObject    handle to check_event_compl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_event_compl


% --- Executes on button press in check_event_buttond.
function check_event_buttond_Callback(hObject, eventdata, handles)
% hObject    handle to check_event_buttond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_event_buttond


% --- Executes on selection change in popup_baseevent.
function popup_baseevent_Callback(hObject, eventdata, handles)
% hObject    handle to popup_baseevent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_baseevent contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_baseevent


% --- Executes during object creation, after setting all properties.
function popup_baseevent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_baseevent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_t_win_stimon_Callback(hObject, eventdata, handles)
% hObject    handle to edit_t_win_stimon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_t_win_stimon as text
%        str2double(get(hObject,'String')) returns contents of edit_t_win_stimon as a double


% --- Executes during object creation, after setting all properties.
function edit_t_win_stimon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_t_win_stimon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_t_win_epoch_Callback(hObject, eventdata, handles)
% hObject    handle to edit_t_win_epoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_t_win_saccade as text
%        str2double(get(hObject,'String')) returns contents of edit_t_win_saccade as a double


% --- Executes during object creation, after setting all properties.
function edit_t_win_epoch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_t_win_epoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_incorrect.
function check_incorrect_Callback(hObject, eventdata, handles)
% hObject    handle to check_incorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_incorrect



function edit_saccade_range_Callback(hObject, eventdata, handles)
% hObject    handle to edit_saccade_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_saccade_range as text
%        str2double(get(hObject,'String')) returns contents of edit_saccade_range as a double


% --- Executes during object creation, after setting all properties.
function edit_saccade_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_saccade_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_set2.
function popup_set2_Callback(hObject, eventdata, handles)
% hObject    handle to popup_set2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_set2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_set2


% --- Executes during object creation, after setting all properties.
function popup_set2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_set2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_event_trialstartd.
function check_event_trialstartd_Callback(hObject, eventdata, handles)
% hObject    handle to check_event_trialstartd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_event_trialstartd


% --- Executes on button press in check_event_trialstartt.
function check_event_trialstartt_Callback(hObject, eventdata, handles)
% hObject    handle to check_event_trialstartt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_event_trialstartt


% --- Executes on button press in check_event_dist1t.
function check_event_dist1t_Callback(hObject, eventdata, handles)
% hObject    handle to check_event_dist1t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_event_dist1t


% --- Executes on button press in check_event_buttont.
function check_event_buttont_Callback(hObject, eventdata, handles)
% hObject    handle to check_event_buttont (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_event_buttont



function edit_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subject as text
%        str2double(get(hObject,'String')) returns contents of edit_subject as a double


% --- Executes during object creation, after setting all properties.
function edit_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_version_Callback(hObject, eventdata, handles)
% hObject    handle to edit_version (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_version as text
%        str2double(get(hObject,'String')) returns contents of edit_version as a double


% --- Executes during object creation, after setting all properties.
function edit_version_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_version (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
