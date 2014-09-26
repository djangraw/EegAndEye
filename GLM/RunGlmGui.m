function varargout = RunGlmGui(varargin)
% RUNGLMGUI M-file for RunGlmGui.fig
%      RUNGLMGUI, by itself, creates a new RUNGLMGUI or raises the existing
%      singleton*.
%
%      H = RUNGLMGUI returns the handle to a new RUNGLMGUI or the handle to
%      the existing singleton*.
%
%      RUNGLMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RUNGLMGUI.M with the given input arguments.
%
%      RUNGLMGUI('Property','Value',...) creates a new RUNGLMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RunGlmGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RunGlmGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Created 7/12 by DJ based on SetUpGlm and RunEegGlm
% Updated 7/30/12 by DJ - sped up by adding events all at once
% Updated 4/30/13 by DJ - responseFns in cells
% Updated 3/19/14 by DJ - tResponse in cells

% Edit the above text to modify the response to help RunGlmGui

% Last Modified by GUIDE v2.5 25-Mar-2013 12:02:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RunGlmGui_OpeningFcn, ...
                   'gui_OutputFcn',  @RunGlmGui_OutputFcn, ...
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


% --- Executes just before RunGlmGui is made visible.
function RunGlmGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RunGlmGui (see VARARGIN)

% Choose default command line output for RunGlmGui
handles.output = hObject;
handles.prefix = 'sq';
handles.experiment = 'Squares';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RunGlmGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RunGlmGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% RUN BUTTON

% --- Executes on button press in button_run.
function button_run_Callback(hObject, eventdata, handles)
% hObject    handle to button_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Clear status
status = {[datestr(now,16) ' - Loading gui parameters...']};
set(handles.list_status,'String',status);
drawnow; % Update view

% Get regressor events
list_events = get(handles.list_events,'String');
nLevels = length(list_events);
regressor_events = cell(nLevels,1);
for iLevel=1:nLevels
    regressor_events{iLevel} = eval(list_events{iLevel});
end
% Get filenames
prefix = get(handles.text_resultsprefix,'String');
suffixes = get(handles.list_filenames, 'String');
filenames = cell(size(suffixes));
for i=1:numel(suffixes)
    if isempty(suffixes{i})
        filenames{i} = '';
    else
        filenames{i} = [prefix, suffixes{i}];
    end
end

% Get settings
datasetprefix = get(handles.text_datasetprefix,'String');
datasetsuffix = get(handles.edit_dataset,'String');
dataset = [datasetprefix datasetsuffix];
offset = str2double(get(handles.edit_offset,'String'));
influence = str2double(get(handles.edit_influence,'String'));
artifact_influence = str2double(get(handles.edit_artifact_influence,'String'));
stddev = str2double(get(handles.edit_stddev,'String'));
vthresh = str2double(get(handles.edit_vthresh,'String'));
all_methods = get(handles.popup_method,'String');
iMethod = get(handles.popup_method,'Value');
method = all_methods{iMethod};

% Get trial rejection rules
trial_rej_rules = {};
trialrej_checkboxes = get(handles.panel_trialrej,'Children'); % gets all the checkboxes in this panel
for i=1:numel(trialrej_checkboxes)
    if get(trialrej_checkboxes(i),'Value')
        trial_rej_rules = [trial_rej_rules, {get(trialrej_checkboxes(i),'String')}];
    end
end

% Get artifacts to reject
artifact_events = {};
artifact_checkboxes = get(handles.panel_artifactrej,'Children'); % gets all the checkboxes in this panel
for i=1:numel(artifact_checkboxes)
    if get(artifact_checkboxes(i),'Value')
        artifact_events = [artifact_events, {get(artifact_checkboxes(i),'String')}];
    end
end


% Update status
status = [status, {sprintf('%s - Preparing dataset %s...',datestr(now,16),dataset)}];
set(handles.list_status,'String',status);
drawnow; % Update view

% Load data
eeglab nogui
if strcmp(dataset(end-3:end),'.set')
    EEG = pop_loadset(dataset);
else
    results = load(dataset);
%     EEG = SubtractOutGlmResponses(results.EEG,results.responseFns,results.regressor_events{results.iLevel});
    EEG = results.EEG;
end

% Add specified events as in SetUpGLM
% load behavioral data
subject = str2double(EEG.subject);
y = loadBehaviorData(subject,[],handles.prefix);
% add events
all_events = [artifact_events regressor_events{:}];
% for iLevel=1:nLevels
%     all_events = [all_events, regressor_events{iLevel}];
% end

% find events
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
EEG = AddEeglabEvents_MultiSession(EEG,y,all_times,all_codes);

% Once all the events are added, add the weights
if ~isfield(EEG.etc,'ureventweights')
    EEG.etc.ureventweights = ones(size(EEG.urevent));
end
for i=1:numel(all_events)
    if ~isempty(weights{i})
        EEG.etc.ureventweights(strcmp(all_events{i},{EEG.urevent(:).type})) = cat(1,weights{i,:});
    end
end


% Epoch data
EEG = pop_epoch(EEG,{'TrialStart-T', 'TrialStart-D'},[-1.5 4.5]);

% Interpolate noisy electrodes and reject noisy trials
EEG = EnforceVoltageThreshold(EEG,vthresh,[regressor_events{:}], artifact_events, influence, artifact_influence);

% Reject bad trials
EEG = RejectEegData(EEG,y,trial_rej_rules);

EEGnew = EEG; % create copy to be modified during run (but only original will be saved)
% Run analysis
for iLevel=1:nLevels % for each level of analysis
    % Update status
    status = [status, {sprintf('%s - Running Analysis Level %d/%d...',datestr(now,16),iLevel,nLevels)}];
    set(handles.list_status,'String',status);
    drawnow; % Update view

    [responseFns{iLevel}, tResponse{iLevel}] = RunEegGlm(EEGnew,regressor_events{iLevel},influence,artifact_events,artifact_influence,method,stddev);
%     responseFns = [responseFns; {newRF}];
%     responseFns = []; tResponse  = [];

    
    if ~isempty(filenames{iLevel})
        % Update status
        status = [status, {sprintf('%s Saving Results as %s...',datestr(now,16), filenames{iLevel})}];
        set(handles.list_status,'String',status);
        drawnow; % Update view
        % Save results
        save(filenames{iLevel},'EEG','responseFns','tResponse','regressor_events',...
            'filenames','iLevel','artifact_events','dataset','offset','influence',...
            'artifact_influence','stddev','vthresh','method','trial_rej_rules');
    end
    EEGnew = SubtractOutGlmResponses(EEGnew,responseFns{iLevel},influence,regressor_events{iLevel},stddev);
end

% Update status
status = [status, {sprintf('%s - Analysis complete!',datestr(now,16))}];
set(handles.list_status,'String',status);
drawnow; % Update view


%% LOAD BUTTON

% --- Executes on button press in button_load.
function button_load_Callback(hObject, eventdata, handles)
% hObject    handle to button_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Pick results file
filename = uigetfile('.mat','Choose GLM Results');
if ~ischar(filename)
    return;
end
results = load(filename);

% Parse dataset
dataset_full = results.dataset; % e.g. 'sq-18-all-filtered-50Hz-noduds.set'
dashes = strfind(dataset_full,'-');
subject_str = dataset_full((dashes(1)+1):(dashes(2)-1));
dataset_suffix = dataset_full((dashes(3)+1):end); 
% Fill in edit boxes
set(handles.edit_subject,'String',subject_str);
set(handles.text_resultsprefix,'String',sprintf('%s-%s-GLMresults-',handles.prefix,subject_str));
set(handles.text_datasetprefix,'String',sprintf('%s-%s-all-',handles.prefix,subject_str));
set(handles.edit_dataset,'String',dataset_suffix);
set(handles.edit_offset,'String',num2str(results.offset));
set(handles.edit_influence,'String',num2str(results.influence));
try set(handles.edit_artifact_influence,'String',num2str(results.artifact_influence));
catch, set(handles.edit_artifact_influence,'String',num2str(results.influence)); % assume influence was same as for events
end  
set(handles.edit_stddev,'String',num2str(results.stddev));
try set(handles.edit_vthresh,'String',num2str(results.vthresh)); 
catch, set(handles.edit_vthresh,'String','???');
end
all_methods = get(handles.popup_method,'String');

% Fill in popups
iMethod = find(strcmp(results.method,all_methods),1);
set(handles.popup_method,'Value',iMethod);

% Fill in checkboxes
% Get trial rejection rules
trialrej_checkboxes = get(handles.panel_trialrej,'Children'); % gets all the checkboxes in this panel
for i=1:numel(trialrej_checkboxes)
    rule = get(trialrej_checkboxes(i),'String');
    isOn = ismember(rule,results.trial_rej_rules);
    set(trialrej_checkboxes(i),'Value',isOn);        
end
% Get artifacts to reject
artifact_checkboxes = get(handles.panel_artifactrej,'Children'); % gets all the checkboxes in this panel
for i=1:numel(artifact_checkboxes)
    event = get(artifact_checkboxes(i),'String');
    isOn = ismember(event,results.artifact_events);
    set(artifact_checkboxes(i),'Value',isOn);        
end

% Fill in Lists
nLevels = length(results.regressor_events);
% Parse events
list_events = cell(nLevels,1);
for iLevel=1:nLevels
    events = results.regressor_events{iLevel};
    eventstr = '{';
    for iEvent = 1:length(events)
        eventstr = [eventstr '''' events{iEvent} ''','];
    end
    eventstr(end) = '}';
    list_events{iLevel} = eventstr;
end
set(handles.list_events,'String',list_events);
% Parse filenames
list_filenames = cell(nLevels,1);
for iLevel = 1:nLevels
    filename = results.filenames{iLevel};
    if isempty(filename)
        list_filenames{iLevel} = '';
    else
        dashes = strfind(filename,'-');    
        list_filenames{iLevel} = filename((dashes(3)+1):end);
    end
end
set(handles.list_filenames,'String',list_filenames);


% Update handles structure
guidata(hObject, handles);

%% OTHER BUTTONS

% --- Executes on button press in button_browse.
function button_browse_Callback(hObject, eventdata, handles)
% hObject    handle to button_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Pick new dataset
filename = uigetfile('*.set;*.mat','Choose starting dataset');
set(handles.edit_dataset,'String',filename);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_add.
function button_add_Callback(hObject, eventdata, handles)
% hObject    handle to button_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get current values
events_list = get(handles.popup_rule,'String');
iEvents = get(handles.popup_rule,'Value');
events = events_list{iEvents};
saveit = get(handles.check_save,'Value');
filename = get(handles.edit_filename,'String');

if ~isempty(events) && iEvents>1
    % Add events
    all_events = get(handles.list_events,'String');
    all_events = [all_events; {events}];
    set(handles.list_events,'String',all_events);
    % Add filename
    all_filenames = get(handles.list_filenames,'String');
    if saveit
        all_filenames = [all_filenames; {filename}];
    else
        all_filenames = [all_filenames; {''}];
    end
    set(handles.list_filenames,'String',all_filenames);
end
   
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_remove.
function button_remove_Callback(hObject, eventdata, handles)
% hObject    handle to button_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get current values
all_events = get(handles.list_events,'String');
all_filenames = get(handles.list_filenames,'String');
iEvents = get(handles.list_events,'Value');

% Remove selected items
all_events(iEvents) = [];
all_filenames(iEvents) = [];
if iEvents>1
    iEvents = iEvents-1;
end

% Update values
set(handles.list_events,'String',all_events);
set(handles.list_events,'Value',iEvents);
set(handles.list_filenames,'String',all_filenames);
set(handles.list_filenames,'Value',iEvents);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_up.
function button_up_Callback(hObject, eventdata, handles)
% hObject    handle to button_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get current values
all_events = get(handles.list_events,'String');
all_filenames = get(handles.list_filenames,'String');
iEvents = get(handles.list_events,'Value');
% Reorder
if iEvents>1
    all_events = all_events([1:iEvents-2, iEvents, iEvents-1, iEvents+1:end]);
    all_filenames = all_filenames([1:iEvents-2, iEvents, iEvents-1, iEvents+1:end]);
    iEvents = iEvents-1;
end
% Update values
set(handles.list_events,'String',all_events);
set(handles.list_events,'Value',iEvents);
set(handles.list_filenames,'String',all_filenames);
set(handles.list_filenames,'Value',iEvents);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in button_down.
function button_down_Callback(hObject, eventdata, handles)
% hObject    handle to button_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get current values
all_events = get(handles.list_events,'String');
all_filenames = get(handles.list_filenames,'String');
iEvents = get(handles.list_events,'Value');
nEvents = length(all_events);
% Reorder
if iEvents<nEvents
    all_events = all_events([1:iEvents-1, iEvents+1, iEvents, iEvents+2:end]);
    all_filenames = all_filenames([1:iEvents-1, iEvents+1, iEvents, iEvents+2:end]);
    iEvents = iEvents+1;
end
% Update values
set(handles.list_events,'String',all_events);
set(handles.list_events,'Value',iEvents);
set(handles.list_filenames,'String',all_filenames);
set(handles.list_filenames,'Value',iEvents);

% Update handles structure
guidata(hObject, handles);



%% LIST BOXES

% --- Executes on selection change in list_events.
function list_events_Callback(hObject, eventdata, handles)
% hObject    handle to list_events (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_events contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_events

% Make sure values for both listboxes match
iEvents = get(hObject,'Value');
set(handles.list_events,'Value',iEvents);
set(handles.list_filenames,'Value',iEvents);
% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in list_filenames.
function list_filenames_Callback(hObject, eventdata, handles)
% hObject    handle to list_filenames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_filenames contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_filenames

% Make sure values for both listboxes match
iEvents = get(hObject,'Value');
set(handles.list_events,'Value',iEvents);
set(handles.list_filenames,'Value',iEvents);
% Update handles structure
guidata(hObject, handles);


%% OTHER

function edit_subject_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subject as text
%        str2double(get(hObject,'String')) returns contents of edit_subject as a double

subject_str = get(hObject,'String');

set(handles.text_resultsprefix,'string',[handles.prefix '-' subject_str '-GLMresults-']);
set(handles.text_datasetprefix,'string',[handles.prefix '-' subject_str '-all-']);
% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in popup_experiment.
function popup_experiment_Callback(hObject, eventdata, handles)
% hObject    handle to popup_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_experiment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_experiment

experiment_cell = get(hObject,'String');
handles.experiment = experiment_cell{get(hObject,'Value')};
switch handles.experiment
    case 'Squares'
        handles.prefix = 'sq';
    case 'SquaresFix'
        handles.prefix = 'sf';
    otherwise
        handles.prefix = '??';      
end

guidata(hObject,handles); % Save variables

% Get subject number
subject_str = get(handles.edit_subject,'string');
% Set prefix for files in Cleanup Panel
set(handles.text_resultsprefix,'string',[handles.prefix '-' subject_str '-GLMresults-']);
set(handles.text_datasetprefix,'string',[handles.prefix '-' subject_str '-all-']);




%% UNUSED

% --- Executes during object creation, after setting all properties.
function list_events_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_events (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function list_filenames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_filenames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_rule.
function popup_rule_Callback(hObject, eventdata, handles)
% hObject    handle to popup_rule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_rule contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_rule


% --- Executes during object creation, after setting all properties.
function popup_rule_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_rule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filename as text
%        str2double(get(hObject,'String')) returns contents of edit_filename as a double


% --- Executes during object creation, after setting all properties.
function edit_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_save.
function check_save_Callback(hObject, eventdata, handles)
% hObject    handle to check_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_save



function edit_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dataset as text
%        str2double(get(hObject,'String')) returns contents of edit_dataset as a double


% --- Executes during object creation, after setting all properties.
function edit_dataset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_offset_Callback(hObject, eventdata, handles)
% hObject    handle to edit_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_offset as text
%        str2double(get(hObject,'String')) returns contents of edit_offset as a double


% --- Executes during object creation, after setting all properties.
function edit_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_influence_Callback(hObject, eventdata, handles)
% hObject    handle to edit_influence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_influence as text
%        str2double(get(hObject,'String')) returns contents of edit_influence as a double


% --- Executes during object creation, after setting all properties.
function edit_influence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_influence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_method.
function popup_method_Callback(hObject, eventdata, handles)
% hObject    handle to popup_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_method


% --- Executes during object creation, after setting all properties.
function popup_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stddev_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stddev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stddev as text
%        str2double(get(hObject,'String')) returns contents of edit_stddev as a double


% --- Executes during object creation, after setting all properties.
function edit_stddev_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stddev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in list_status.
function list_status_Callback(hObject, eventdata, handles)
% hObject    handle to list_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_status contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_status


% --- Executes during object creation, after setting all properties.
function list_status_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_cross.
function check_cross_Callback(hObject, eventdata, handles)
% hObject    handle to check_cross (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_cross


% --- Executes on button press in check_errant.
function check_errant_Callback(hObject, eventdata, handles)
% hObject    handle to check_errant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_errant


% --- Executes on button press in check_blinkstart.
function check_blinkstart_Callback(hObject, eventdata, handles)
% hObject    handle to check_blinkstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_blinkstart


% --- Executes on button press in check_blinkend.
function check_blinkend_Callback(hObject, eventdata, handles)
% hObject    handle to check_blinkend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_blinkend


% --- Executes on button press in check_button.
function check_button_Callback(hObject, eventdata, handles)
% hObject    handle to check_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_button


% --- Executes on button press in check_backward_trial.
function check_backward_trial_Callback(hObject, eventdata, handles)
% hObject    handle to check_backward_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_backward_trial


% --- Executes on button press in check_skip_trial.
function check_skip_trial_Callback(hObject, eventdata, handles)
% hObject    handle to check_skip_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_skip_trial


% --- Executes on button press in check_skippedends_trial.
function check_skippedends_trial_Callback(hObject, eventdata, handles)
% hObject    handle to check_skippedends_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_skippedends_trial


% --- Executes on button press in check_errant_trial.
function check_errant_trial_Callback(hObject, eventdata, handles)
% hObject    handle to check_errant_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_errant_trial


% --- Executes on button press in check_blink_trial.
function check_blink_trial_Callback(hObject, eventdata, handles)
% hObject    handle to check_blink_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_blink_trial


% --- Executes on button press in check_earlybutton_trial.
function check_earlybutton_trial_Callback(hObject, eventdata, handles)
% hObject    handle to check_earlybutton_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_earlybutton_trial


% --- Executes on button press in check_latebutton_trial.
function check_latebutton_trial_Callback(hObject, eventdata, handles)
% hObject    handle to check_latebutton_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_latebutton_trial



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



function edit_vthresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vthresh as text
%        str2double(get(hObject,'String')) returns contents of edit_vthresh as a double


% --- Executes during object creation, after setting all properties.
function edit_vthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popup_experiment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_artifact_influence_Callback(hObject, eventdata, handles)
% hObject    handle to edit_artifact_influence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_artifact_influence as text
%        str2double(get(hObject,'String')) returns contents of edit_artifact_influence as a double


% --- Executes during object creation, after setting all properties.
function edit_artifact_influence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_artifact_influence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_wrongbutton_trial.
function check_wrongbutton_trial_Callback(hObject, eventdata, handles)
% hObject    handle to check_wrongbutton_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_wrongbutton_trial
