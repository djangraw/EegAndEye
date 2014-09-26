function varargout = GetLooResults(varargin)
% GETLOORESULTS M-file for GetLooResults.fig
%      GETLOORESULTS, by itself, creates a new GETLOORESULTS or raises the existing
%      singleton*.
%
%      H = GETLOORESULTS returns the handle to a new GETLOORESULTS or the handle to
%      the existing singleton*.
%
%      GETLOORESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GETLOORESULTS.M with the given input arguments.
%
%      GETLOORESULTS('Property','Value',...) creates a new GETLOORESULTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GetLooResults_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GetLooResults_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GetLooResults

% Last Modified by GUIDE v2.5 24-Feb-2011 16:27:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GetLooResults_OpeningFcn, ...
                   'gui_OutputFcn',  @GetLooResults_OutputFcn, ...
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


% --- Executes just before GetLooResults is made visible.
function GetLooResults_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GetLooResults (see VARARGIN)

% Choose default command line output for GetLooResults
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GetLooResults wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GetLooResults_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     IMPORTANT FUNCTIONS                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in button_load.
function button_load_Callback(hObject, eventdata, handles)
% hObject    handle to button_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.filename = get(handles.edit_filename,'String');
if ~exist(handles.filename)
    set(handles.text_file,'String','File Not Found!');
else
    load(handles.filename)
    set(handles.text_file,'String',sprintf('%d entries found!',length(LOO)));
    set(handles.popup_setname1,'String',[{'Any'}, unique({LOO(:).setname1}) ] );
    set(handles.popup_setname2,'String',[{'Any'}, unique({LOO(:).setname2}) ] );
    set(handles.popup_reference,'String',[{'Any'}, num2str(unique(cellfun('length',{LOO(:).reference}))') ] );
    set(handles.popup_twl,'String',[{'Any'}, num2str(unique([LOO(:).trainingwindowlength])') ] );
    set(handles.popup_twi,'String',[{'Any'}, num2str(unique([LOO(:).trainingwindowinterval])') ] );
    % etcetera...
end
guidata(hObject, handles);

% --- Executes on button press in button_search.
function button_search_Callback(hObject, eventdata, handles)
% hObject    handle to button_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Parse input from ChooseSearchCriteria panel of GUI
inputStr = ['''filename'', ''' handles.filename ''''];
dateStr = get(handles.edit_date,'String');
if ~isempty(dateStr)
    inputStr = [inputStr, ', ''datetime'', ''' dateStr ''''];
end
setname1Str = get(handles.popup_setname1,'String');
setname1Str = setname1Str{get(handles.popup_setname1,'Value')};
if ~strcmp(setname1Str,'Any')
    inputStr = [inputStr ', ''setname1'', ''' setname1Str ''''];
end
setname2Str = get(handles.popup_setname2,'String');
setname2Str = setname2Str{get(handles.popup_setname2,'Value')};
if ~strcmp(setname2Str,'Any')
    inputStr = [inputStr ', ''setname2'', ''' setname2Str ''''];
end
referenceStr = get(handles.popup_reference,'String');
referenceStr = referenceStr{get(handles.popup_reference,'Value')};
if ~strcmp(referenceStr,'Any')
    inputStr = [inputStr ', ''reference'', ''' referenceStr ''''];
end
twlStr = get(handles.popup_twl,'String');
twlStr = twlStr{get(handles.popup_twl,'Value')};
if ~strcmp(twlStr,'Any')
    inputStr = [inputStr ', ''trainingwindowlength'', ''' twlStr ''''];
end
twiStr = get(handles.popup_twi,'String');
twiStr = twiStr{get(handles.popup_twi,'Value')};
if ~strcmp(twiStr,'Any')
    inputStr = [inputStr ', ''trainingwindowinterval'', ''' twiStr ''''];
end

% Get results from LOO struct
eval(['handles.fits = LoadLooResults(' inputStr ');']);
% Send to PlotResults panel of GUI
if isempty(handles.fits)
    set(handles.popup_results,'String',' ');
else
    set(handles.popup_results,'String',{handles.fits.datetime});
end
set(handles.popup_results,'Value',1);
set(handles.text_results,'String',sprintf('%d result(s) found: pick one above',numel(handles.fits)));
guidata(hObject, handles);


% --- Executes on button press in button_calendar.
function button_calendar_Callback(hObject, eventdata, handles)
% hObject    handle to button_calendar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uicalendar('Weekend', [1 0 0 0 0 0 1], ...  
'SelectionType', 1, ...  
'DestinationUI', handles.edit_date);



% --- Executes on selection change in popup_results.
function popup_results_Callback(hObject, eventdata, handles)
% hObject    handle to popup_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_results contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        popup_results

if isfield(handles,'fits') && ~isempty(handles.fits)
    result = handles.fits(get(handles.popup_results,'Value'));
    fields = fieldnames(result);
    dispText = {''};
    for i=1:numel(fields)
        thisField = result.(fields{i});
        if isnumeric(thisField) || islogical(thisField)
            if length(thisField)==1
                dispText = [dispText; {sprintf('%s: %g',fields{i},thisField)}];
            else
                dispText = [dispText; {sprintf('%s: %g to %g',fields{i},thisField(1),thisField(end))}];
            end
        else
            dispText = [dispText; {sprintf('%s: %s',fields{i},thisField)}];
        end
    end
    set(handles.text_results,'String',dispText);
    guidata(hObject, handles);
end


% --- Executes on button press in button_plot.
function button_plot_Callback(hObject, eventdata, handles)
% hObject    handle to button_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'fits') && ~isempty(handles.fits)
    result = handles.fits(get(handles.popup_results,'Value'));
    figure;
    PlotLooResults(result);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     OTHER FUNCTIONS                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


function edit_date_Callback(hObject, eventdata, handles)
% hObject    handle to edit_date (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_date as text
%        str2double(get(hObject,'String')) returns contents of edit_date as a double


% --- Executes during object creation, after setting all properties.
function edit_date_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_date (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_setname1.
function popup_setname1_Callback(hObject, eventdata, handles)
% hObject    handle to popup_setname1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_setname1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_setname1


% --- Executes during object creation, after setting all properties.
function popup_setname1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_setname1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_setname2.
function popup_setname2_Callback(hObject, eventdata, handles)
% hObject    handle to popup_setname2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_setname2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_setname2


% --- Executes during object creation, after setting all properties.
function popup_setname2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_setname2 (see GCBO)
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


% --- Executes on selection change in popup_twl.
function popup_twl_Callback(hObject, eventdata, handles)
% hObject    handle to popup_twl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_twl contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_twl


% --- Executes during object creation, after setting all properties.
function popup_twl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_twl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_twi.
function popup_twi_Callback(hObject, eventdata, handles)
% hObject    handle to popup_twi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_twi contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_twi


% --- Executes during object creation, after setting all properties.
function popup_twi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_twi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popup_results_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
