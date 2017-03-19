function varargout = SkyrmionsControlPanel(varargin)
% SKYRMIONSCONTROLPANEL M-file for SkyrmionsControlPanel.fig
%      SKYRMIONSCONTROLPANEL, by itself, creates a new SKYRMIONSCONTROLPANEL or raises the existing
%      singleton*.
%
%      H = SKYRMIONSCONTROLPANEL returns the handle to a new SKYRMIONSCONTROLPANEL or the handle to
%      the existing singleton*.
%
%      SKYRMIONSCONTROLPANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SKYRMIONSCONTROLPANEL.M with the given input arguments.
%
%      SKYRMIONSCONTROLPANEL('Property','Value',...) creates a new SKYRMIONSCONTROLPANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SkyrmionsControlPanel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SkyrmionsControlPanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SkyrmionsControlPanel

% Last Modified by GUIDE v2.5 01-Aug-2012 12:27:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SkyrmionsControlPanel_OpeningFcn, ...
                   'gui_OutputFcn',  @SkyrmionsControlPanel_OutputFcn, ...
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


% --- Executes just before SkyrmionsControlPanel is made visible.
function SkyrmionsControlPanel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SkyrmionsControlPanel (see VARARGIN)

% Choose default command line output for SkyrmionsControlPanel
handles.output = hObject;

set(handles.PlottingOptionsButtonGroup,'SelectionChangeFcn',@PlottingOptionsButtonGroup_SelectionChangeFcn);
set(handles.WhatToPlotButtonGroup,'SelectionChangeFcn',@WhatToPlotButtonGroup_SelectionChangeFcn);
set(handles.IntegralFixedAxisButtonGroup,'SelectionChangeFcn',@IntegralFixedAxisButtonGroup_SelectionChangeFcn);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SkyrmionsControlPanel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SkyrmionsControlPanel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton2



function SystemParametersFileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to SystemParametersFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SystemParametersFileEdit as text
%        str2double(get(hObject,'String')) returns contents of SystemParametersFileEdit as a double
filenameParameters=get(hObject,'String');
if (isempty(filenameParameters))
    set(hObject,'String','Enter path to the file containing the system parameters');
    errordlg('Could not read filename.','System Parameters');
else
    fileID=fopen(regexprep(filenameParameters,' ',''));
    if (fileID==-1)
        errordlg(['Could not open file ',filenameParameters],'System Parameters');
    else
        fclose(fileID);
        filenameSysParams=setSysParamsFilename(filenameParameters);
        getSysParams(filenameParameters,filenameSysParams);
        load(filenameSysParams);
        
        % Setting ImaginaryTimesPopUpMenu entries
        jmax=size(SysParams__saving_times_in_filenames);
        jmax=jmax(1);
        ImaginaryTimesPopUpMenu_string=SysParams__saving_times_in_filenames(1,:);
        for j=2:jmax
            ImaginaryTimesPopUpMenu_string=[ImaginaryTimesPopUpMenu_string,'|',SysParams__saving_times_in_filenames(j,:)];
        end
        h_ImaginaryTimesPopUpMenu=findobj('Tag','ImaginaryTimesPopUpMenu');
        set(h_ImaginaryTimesPopUpMenu,'String',ImaginaryTimesPopUpMenu_string);
        clear h_ImaginaryTimesPopUpMenu ImaginaryTimesPopUpMenu_string ...
              jmax j;
        
        % Setting plotting options
        % Fixed x, Fixed y, Fixed z
        h=findobj('Tag','FixedXMinXText');
        set(h,'String',num2str(-0.5*(SysParams__Mx-1)*SysParams__ax,'%3.4f'));
        h=findobj('Tag','FixedXMaxXText');
        set(h,'String',num2str( 0.5*(SysParams__Mx-1)*SysParams__ax,'%3.4f'));
        h=findobj('Tag','FixedYMinYText');
        set(h,'String',num2str(-0.5*(SysParams__My-1)*SysParams__ay,'%3.4f'));
        h=findobj('Tag','FixedYMaxYText');
        set(h,'String',num2str( 0.5*(SysParams__My-1)*SysParams__ay,'%3.4f'));
        h=findobj('Tag','FixedZMinZText');
        set(h,'String',num2str(-0.5*(SysParams__Mz-1)*SysParams__az,'%3.4f'));
        h=findobj('Tag','FixedZMaxZText');
        set(h,'String',num2str( 0.5*(SysParams__Mz-1)*SysParams__az,'%3.4f'));
        
        h=findobj('Tag','FixedXEdit');
        set(h,'String',num2str(0));
        h=findobj('Tag','FixedYEdit');
        set(h,'String',num2str(0));
        h=findobj('Tag','FixedZEdit');
        set(h,'String',num2str(0));
        
        h=findobj('Tag','FixedXSlider');
        set(h,'Min',-0.5*(SysParams__Mx-1)*SysParams__ax);
        set(h,'Max', 0.5*(SysParams__Mx-1)*SysParams__ax);
        h=findobj('Tag','FixedYSlider');
        set(h,'Min',-0.5*(SysParams__My-1)*SysParams__ay);
        set(h,'Max', 0.5*(SysParams__My-1)*SysParams__ay);
        h=findobj('Tag','FixedZSlider');
        set(h,'Min',-0.5*(SysParams__Mz-1)*SysParams__az);
        set(h,'Max', 0.5*(SysParams__Mz-1)*SysParams__az);
        
        % Ellipsoid
        h=findobj('Tag','EllipsoidMaxRxText');
        set(h,'String',num2str(0.5*(SysParams__Mx-1)*SysParams__ax,'%3.4f'));
        h=findobj('Tag','EllipsoidMaxRyText');
        set(h,'String',num2str(0.5*(SysParams__My-1)*SysParams__ay,'%3.4f'));
        h=findobj('Tag','EllipsoidMaxRzText');
        set(h,'String',num2str(0.5*(SysParams__Mz-1)*SysParams__az,'%3.4f'));
        
        h=findobj('Tag','EllipsoidRxSlider');
        minValue=0;
        maxValue=0.5*(SysParams__Mx-1)*SysParams__ax;
        sliderValue=(minValue+maxValue)*0.5;
        set(h,'Value',sliderValue);
        set(h,'Min',minValue);
        set(h,'Max',maxValue);
        h=findobj('Tag','EllipsoidRxEdit');
        set(h,'String',num2str(sliderValue,'%3.4f'));
        
        h=findobj('Tag','EllipsoidRySlider');
        minValue=0;
        maxValue=0.5*(SysParams__My-1)*SysParams__ay;
        sliderValue=(minValue+maxValue)*0.5;
        set(h,'Value',sliderValue);
        set(h,'Min',minValue);
        set(h,'Max',maxValue);
        h=findobj('Tag','EllipsoidRyEdit');
        set(h,'String',num2str(sliderValue,'%3.4f'));

        h=findobj('Tag','EllipsoidRzSlider');
        minValue=0;
        maxValue=0.5*(SysParams__Mz-1)*SysParams__az;
        sliderValue=(minValue+maxValue)*0.5;
        set(h,'Value',sliderValue);
        set(h,'Min',minValue);
        set(h,'Max',maxValue);
        h=findobj('Tag','EllipsoidRzEdit');
        set(h,'String',num2str(sliderValue,'%3.4f'));    
        
        % Line
        set(handles.LineWidthSlider,'Max',(SysParams__Mz-1)*SysParams__az);
        set(handles.LineWidthSlider,'Value',(SysParams__Mz-1)*SysParams__az);
        set(handles.LineWidthMaxWidthText,'String',num2str((SysParams__Mz-1)*SysParams__az,'%3.3f'));
        set(handles.LineWidthEdit,'String',num2str((SysParams__Mz-1)*SysParams__az,'%3.3f'));
        set(handles.LineWidthEdit,'String',num2str((SysParams__Mz-1)*SysParams__az,'%3.3f'));
    end
end




% --- Executes during object creation, after setting all properties.
function SystemParametersFileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SystemParametersFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ImaginaryTimesPopUpMenu.
function ImaginaryTimesPopUpMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ImaginaryTimesPopUpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ImaginaryTimesPopUpMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImaginaryTimesPopUpMenu

if(get(handles.IsosurfaceRadioButton,'Value'))
        filenameParameters=get(handles.SystemParametersFileEdit,'String');
        filenameSysParams=setSysParamsFilename(filenameParameters);
        load(filenameSysParams);
        j=get(handles.ImaginaryTimesPopUpMenu,'Value');
        filenamePsis=SysParams__filenames(j,:);
        % Get a low resolution sample of the data, to estimate the minimal 
        % and maximal values.
        intersp_x=min([SysParams__Mx,max([1,round(SysParams__Mx/64)])]);
        intersp_y=min([SysParams__My,max([1,round(SysParams__My/64)])]);
        intersp_z=min([SysParams__Mz,max([1,round(SysParams__Mz/64)])]);
        [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,x,y,z]=...
         getPsisInterspaced(filenamePsis,filenameSysParams,...
                            intersp_x,intersp_y,intersp_z,...
                            get(handles.hedgehogizeCheckBox,'Value')*(1+get(handles.DiffFromHedgehogCheckebox,'Value')));
        % Chose which datatype is selected in WhatToPlotButtonGroup
        if(get(handles.AbsPsiRadioButton,'Value'))
            WhatToPlot=get(handles.AbsPsiRadioButton,'String');
        elseif(get(handles.ArgPsiRadioButton,'Value'))
            WhatToPlot=get(handles.ArgPsiRadioButton,'String');
        elseif(get(handles.RePsiRadioButton,'Value'))
            WhatToPlot=get(handles.RePsiRadioButton,'String');
        elseif(get(handles.ImPsiRadioButton,'Value'))
            WhatToPlot=get(handles.ImPsiRadioButton,'String');
        elseif(get(handles.nRadioButton,'Value'))
            WhatToPlot=get(handles.nRadioButton,'String');
        elseif(get(handles.nTOFRadioButton,'Value'))
            WhatToPlot=get(handles.nTOFRadioButton,'String');
        elseif(get(handles.AbsPsi_2RadioButton,'Value'))
            WhatToPlot=get(handles.AbsPsi_2RadioButton,'String');
        end
        WhatToPlot=strtrim(WhatToPlot);
        ComponentNo=get(handles.WhatToPlotPopUpMenu,'Value');
        Component=get(handles.WhatToPlotPopUpMenu,'String');
        Component=char(Component(ComponentNo));
        
        if get(handles.HedgehogBasisRadioButton,'Value')==1
            whichBasis=1;
        else
            whichBasis=2;
        end
        ArrayToPlot=setArrayToPlot(Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,...
                                   nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,...
                                   WhatToPlot,Component,whichBasis);
        s=size(ArrayToPlot);
        jmax=size(s); jmax=jmax(2);
        minValue=min(ArrayToPlot);
        maxValue=max(ArrayToPlot);
        for j=2:jmax
            minValue=min(minValue);
            maxValue=max(maxValue);
        end
        if(minValue==maxValue)
            minValue=minValue-1e-5*abs(minValue);
            maxValue=maxValue+1e-5*abs(maxValue);
            if(minValue==0)
                minValue=-1e-5;
                maxValue=1e-5;
            end
        end
        h=findobj('Tag','IsosurfaceMinText');
        set(h,'String',num2str(minValue,'%3.4f'));
        h=findobj('Tag','IsosurfaceMaxText');
        set(h,'String',num2str(maxValue,'%3.4f'));
        h=findobj('Tag','IsosurfaceSlider');
        set(h,'Min',minValue);
        set(h,'Max',maxValue);
        set(h,'Value',0.5*(minValue+maxValue));
        set(handles.IsosurfaceEdit,'String',num2str(0.5*(minValue+maxValue)));
    clear  Psi1 Psi2 Psi3 n1 n2 n3 nTOF1 nTOF2 nTOF3 ...
           ArrayToPlot filenameParameters filenameSysParams ...
           j s jmax dataTypeName intersp_x intersp_y intersp_z;
end


% --- Executes during object creation, after setting all properties.
function ImaginaryTimesPopUpMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImaginaryTimesPopUpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function FixedXSlider_Callback(hObject, eventdata, handles)
% hObject    handle to FixedXSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderValue=get(handles.FixedXSlider,'Value');
h=findobj('Tag','FixedXEdit');
set(h,'String',num2str(sliderValue,'%3.4f'));



% --- Executes during object creation, after setting all properties.
function FixedXSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FixedXSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




function FixedXEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FixedXEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FixedXEdit as text
%        str2double(get(hObject,'String')) returns contents of FixedXEdit as a double

editValue=str2num(get(handles.FixedXEdit,'String'));
if (editValue>get(handles.FixedXSlider,'Max'))
    editValue=get(handles.FixedXSlider,'Max');
    set(handles.FixedXEdit,'String',num2str(editValue,'%3.4f'));
elseif (editValue<get(handles.FixedXSlider,'Min'))
    editValue=get(handles.FixedXSlider,'Min');
    set(handles.FixedXEdit,'String',num2str(editValue,'%3.4f'));
end
set(handles.FixedXSlider,'Value',editValue);


% --- Executes during object creation, after setting all properties.
function FixedXEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FixedXEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function FixedYSlider_Callback(hObject, eventdata, handles)
% hObject    handle to FixedYSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


sliderValue=get(handles.FixedYSlider,'Value');
h=findobj('Tag','FixedYEdit');
set(h,'String',num2str(sliderValue,'%3.4f'));


% --- Executes during object creation, after setting all properties.
function FixedYSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FixedYSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function FixedYEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FixedYEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FixedYEdit as text
%        str2double(get(hObject,'String')) returns contents of FixedYEdit as a double

editValue=str2num(get(handles.FixedYEdit,'String'));
if (editValue>get(handles.FixedYSlider,'Max'))
    editValue=get(handles.FixedYSlider,'Max');
    set(handles.FixedYEdit,'String',num2str(editValue,'%3.4f'));
elseif (editValue<get(handles.FixedYSlider,'Min'))
    editValue=get(handles.FixedYSlider,'Min');
    set(handles.FixedYEdit,'String',num2str(editValue,'%3.4f'));
end
set(handles.FixedYSlider,'Value',editValue);


% --- Executes during object creation, after setting all properties.
function FixedYEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FixedYEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function FixedZSlider_Callback(hObject, eventdata, handles)
% hObject    handle to FixedZSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


sliderValue=get(handles.FixedZSlider,'Value');
h=findobj('Tag','FixedZEdit');
set(h,'String',num2str(sliderValue,'%3.4f'));


% --- Executes during object creation, after setting all properties.
function FixedZSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FixedZSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function FixedZEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FixedZEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FixedZEdit as text
%        str2double(get(hObject,'String')) returns contents of FixedZEdit as a double

editValue=str2num(get(handles.FixedZEdit,'String'));
if (editValue>get(handles.FixedZSlider,'Max'))
    editValue=get(handles.FixedZSlider,'Max');
    set(handles.FixedZEdit,'String',num2str(editValue,'%3.4f'));
elseif (editValue<get(handles.FixedZSlider,'Min'))
    editValue=get(handles.FixedZSlider,'Min');
    set(handles.FixedZEdit,'String',num2str(editValue,'%3.4f'));
end
set(handles.FixedZSlider,'Value',editValue);


% --- Executes during object creation, after setting all properties.
function FixedZEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FixedZEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function EllipsoidRxSlider_Callback(hObject, eventdata, handles)
% hObject    handle to EllipsoidRxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderValue=get(handles.EllipsoidRxSlider,'Value');
if (get(handles.FixedRatioCheckBox,'Value'))
    rxPrevValue=str2num(get(handles.EllipsoidRxEdit,'String'));
    ryPrevValue=str2num(get(handles.EllipsoidRyEdit,'String'));
    rzPrevValue=str2num(get(handles.EllipsoidRzEdit,'String'));
    rxMax=get(handles.EllipsoidRxSlider,'Max');
    ryMax=get(handles.EllipsoidRySlider,'Max');
    rzMax=get(handles.EllipsoidRzSlider,'Max');
    if(rxPrevValue~=0)
        mult=sliderValue/rxPrevValue;
        if(mult*rxPrevValue>=rxMax || mult*ryPrevValue>=ryMax || mult*rzPrevValue>=rzMax )
            mult=min([rxMax/rxPrevValue,ryMax/ryPrevValue,rzMax/rzPrevValue]);
        end
        rxNewValue=mult*rxPrevValue;
        ryNewValue=mult*ryPrevValue;
        rzNewValue=mult*rzPrevValue;
    else
        rxNewValue=rxMax*0.0001;
        ryNewValue=ryPrevValue;
        rzNewValue=rzPrevValue;
        if(ryPrevValue==0)
            ryNewValue=ryMax*0.0001;
        end
        if(rzPrevValue==0)
            rzNewValue=rzMax*0.0001;
        end
    end
    set(handles.EllipsoidRxSlider,'Value',rxNewValue);
    set(handles.EllipsoidRySlider,'Value',ryNewValue);
    set(handles.EllipsoidRzSlider,'Value',rzNewValue);
    set(handles.EllipsoidRxEdit,'String',num2str(rxNewValue,'%3.4f'));
    set(handles.EllipsoidRyEdit,'String',num2str(ryNewValue,'%3.4f'));
    set(handles.EllipsoidRzEdit,'String',num2str(rzNewValue,'%3.4f'));
else
    h=findobj('Tag','EllipsoidRxEdit');
    set(h,'String',num2str(sliderValue,'%3.4f'));
end


% --- Executes during object creation, after setting all properties.
function EllipsoidRxSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EllipsoidRxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function EllipsoidRxEdit_Callback(hObject, eventdata, handles)
% hObject    handle to EllipsoidRxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EllipsoidRxEdit as text
%        str2double(get(hObject,'String')) returns contents of EllipsoidRxEdit as a double

editValue=str2num(get(handles.EllipsoidRxEdit,'String'));
if (get(handles.FixedRatioCheckBox,'Value'))
    rxPrevValue=get(handles.EllipsoidRxSlider,'Value');
    ryPrevValue=get(handles.EllipsoidRySlider,'Value');
    rzPrevValue=get(handles.EllipsoidRzSlider,'Value');
    rxMax=get(handles.EllipsoidRxSlider,'Max');
    ryMax=get(handles.EllipsoidRySlider,'Max');
    rzMax=get(handles.EllipsoidRzSlider,'Max');
    rxMin=get(handles.EllipsoidRxSlider,'Min');
    ryMin=get(handles.EllipsoidRySlider,'Min');
    rzMin=get(handles.EllipsoidRzSlider,'Min');
    if(rxPrevValue~=0)
        mult=editValue/rxPrevValue;
        if(mult*rxPrevValue>=rxMax || mult*ryPrevValue>=ryMax || mult*rzPrevValue>=rzMax )
            mult=min([rxMax/rxPrevValue,ryMax/ryPrevValue,rzMax/rzPrevValue]);
        end
        if(mult*rxPrevValue<=rxMin || mult*ryPrevValue<=ryMin || mult*rzPrevValue<=rzMin)
            mult=max([rxMin/rxPrevValue,ryMin/ryPrevValue,rzMin/rzPrevValue]);
        end
        rxNewValue=mult*rxPrevValue;
        ryNewValue=mult*ryPrevValue;
        rzNewValue=mult*rzPrevValue;
    else
        rxNewValue=rxMax*0.0001;
        ryNewValue=ryPrevValue;
        rzNewValue=rzPrevValue;
        if(ryPrevValue==0)
            ryNewValue=ryMax*0.0001;
        end
        if(rzPrevValue==0)
            rzNewValue=rzMax*0.0001;
        end
    end
    set(handles.EllipsoidRxSlider,'Value',rxNewValue);
    set(handles.EllipsoidRySlider,'Value',ryNewValue);
    set(handles.EllipsoidRzSlider,'Value',rzNewValue);
    set(handles.EllipsoidRxEdit,'String',num2str(rxNewValue,'%3.4f'));
    set(handles.EllipsoidRyEdit,'String',num2str(ryNewValue,'%3.4f'));
    set(handles.EllipsoidRzEdit,'String',num2str(rzNewValue,'%3.4f'));
else
    if (editValue>get(handles.EllipsoidRxSlider,'Max'))
        editValue=get(handles.EllipsoidRxSlider,'Max');
        set(handles.EllipsoidRxEdit,'String',num2str(editValue,'%3.4f'));
    elseif (editValue<get(handles.EllipsoidRxSlider,'Min'))
        editValue=get(handles.EllipsoidRxSlider,'Min');
        set(handles.EllipsoidRxEdit,'String',num2str(editValue,'%3.4f'));
    end
    set(handles.EllipsoidRxSlider,'Value',editValue);
end


% --- Executes during object creation, after setting all properties.
function EllipsoidRxEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EllipsoidRxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function EllipsoidRySlider_Callback(hObject, eventdata, handles)
% hObject    handle to EllipsoidRySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderValue=get(handles.EllipsoidRySlider,'Value');
if (get(handles.FixedRatioCheckBox,'Value'))
    rxPrevValue=str2num(get(handles.EllipsoidRxEdit,'String'));
    ryPrevValue=str2num(get(handles.EllipsoidRyEdit,'String'));
    rzPrevValue=str2num(get(handles.EllipsoidRzEdit,'String'));
    rxMax=get(handles.EllipsoidRxSlider,'Max');
    ryMax=get(handles.EllipsoidRySlider,'Max');
    rzMax=get(handles.EllipsoidRzSlider,'Max');
    if(ryPrevValue~=0)
        mult=sliderValue/ryPrevValue;
        if(mult*rxPrevValue>=rxMax || mult*ryPrevValue>=ryMax || mult*rzPrevValue>=rzMax )
            mult=min([rxMax/rxPrevValue,ryMax/ryPrevValue,rzMax/rzPrevValue]);
        end
        rxNewValue=mult*rxPrevValue;
        ryNewValue=mult*ryPrevValue;
        rzNewValue=mult*rzPrevValue;
    else
        ryNewValue=ryMax*0.0001;
        rxNewValue=rxPrevValue;
        rzNewValue=rzPrevValue;
        if(rxPrevValue==0)
            rxNewValue=rxMax*0.0001;
        end
        if(rzPrevValue==0)
            rzNewValue=rzMax*0.0001;
        end
    end
    set(handles.EllipsoidRxSlider,'Value',rxNewValue);
    set(handles.EllipsoidRySlider,'Value',ryNewValue);
    set(handles.EllipsoidRzSlider,'Value',rzNewValue);
    set(handles.EllipsoidRxEdit,'String',num2str(rxNewValue,'%3.4f'));
    set(handles.EllipsoidRyEdit,'String',num2str(ryNewValue,'%3.4f'));
    set(handles.EllipsoidRzEdit,'String',num2str(rzNewValue,'%3.4f'));
else
    h=findobj('Tag','EllipsoidRyEdit');
    set(h,'String',num2str(sliderValue,'%3.4f'));
end


% --- Executes during object creation, after setting all properties.
function EllipsoidRySlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EllipsoidRySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function EllipsoidRyEdit_Callback(hObject, eventdata, handles)
% hObject    handle to EllipsoidRyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EllipsoidRyEdit as text
%        str2double(get(hObject,'String')) returns contents of EllipsoidRyEdit as a double

editValue=str2num(get(handles.EllipsoidRyEdit,'String'));
if (editValue>get(handles.EllipsoidRySlider,'Max'))
    editValue=get(handles.EllipsoidRySlider,'Max');
    set(handles.EllipsoidRyEdit,'String',num2str(editValue,'%3.4f'));
elseif (editValue<get(handles.EllipsoidRySlider,'Min'))
    editValue=get(handles.EllipsoidRySlider,'Min');
    set(handles.EllipsoidRyEdit,'String',num2str(editValue,'%3.4f'));
end
set(handles.EllipsoidRySlider,'Value',editValue);


% --- Executes during object creation, after setting all properties.
function EllipsoidRyEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EllipsoidRyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function EllipsoidRzSlider_Callback(hObject, eventdata, handles)
% hObject    handle to EllipsoidRzSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderValue=get(handles.EllipsoidRzSlider,'Value');
if (get(handles.FixedRatioCheckBox,'Value'))
    rxPrevValue=str2num(get(handles.EllipsoidRxEdit,'String'));
    ryPrevValue=str2num(get(handles.EllipsoidRyEdit,'String'));
    rzPrevValue=str2num(get(handles.EllipsoidRzEdit,'String'));
    rxMax=get(handles.EllipsoidRxSlider,'Max');
    ryMax=get(handles.EllipsoidRySlider,'Max');
    rzMax=get(handles.EllipsoidRzSlider,'Max');
    if(rzPrevValue~=0)
        mult=sliderValue/rzPrevValue;
        if(mult*rxPrevValue>=rxMax || mult*ryPrevValue>=ryMax || mult*rzPrevValue>=rzMax )
            mult=min([rxMax/rxPrevValue,ryMax/ryPrevValue,rzMax/rzPrevValue]);
        end
        rxNewValue=mult*rxPrevValue;
        ryNewValue=mult*ryPrevValue;
        rzNewValue=mult*rzPrevValue;
    else
        rzNewValue=rzMax*0.0001;
        rxNewValue=rxPrevValue;
        ryNewValue=ryPrevValue;
        if(rxPrevValue==0)
            rxNewValue=rxMax*0.0001;
        end
        if(ryPrevValue==0)
            ryNewValue=ryMax*0.0001;
        end
    end
    set(handles.EllipsoidRxSlider,'Value',rxNewValue);
    set(handles.EllipsoidRySlider,'Value',ryNewValue);
    set(handles.EllipsoidRzSlider,'Value',rzNewValue);
    set(handles.EllipsoidRxEdit,'String',num2str(rxNewValue,'%3.4f'));
    set(handles.EllipsoidRyEdit,'String',num2str(ryNewValue,'%3.4f'));
    set(handles.EllipsoidRzEdit,'String',num2str(rzNewValue,'%3.4f'));
else
    h=findobj('Tag','EllipsoidRyEdit');
    set(h,'String',num2str(sliderValue,'%3.4f'));
end


% --- Executes during object creation, after setting all properties.
function EllipsoidRzSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EllipsoidRzSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function EllipsoidRzEdit_Callback(hObject, eventdata, handles)
% hObject    handle to EllipsoidRzEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EllipsoidRzEdit as text
%        str2double(get(hObject,'String')) returns contents of EllipsoidRzEdit as a double

editValue=str2num(get(handles.EllipsoidRzEdit,'String'));
if (editValue>get(handles.EllipsoidRzSlider,'Max'))
    editValue=get(handles.EllipsoidRzSlider,'Max');
    set(handles.EllipsoidRzEdit,'String',num2str(editValue,'%3.4f'));
elseif (editValue<get(handles.EllipsoidRzSlider,'Min'))
    editValue=get(handles.EllipsoidRzSlider,'Min');
    set(handles.EllipsoidRzEdit,'String',num2str(editValue,'%3.4f'));
end
set(handles.EllipsoidRzSlider,'Value',editValue);


% --- Executes during object creation, after setting all properties.
function EllipsoidRzEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EllipsoidRzEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function IsosurfaceSlider_Callback(hObject, eventdata, handles)
% hObject    handle to IsosurfaceSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderValue=get(handles.IsosurfaceSlider,'Value');
h=findobj('Tag','IsosurfaceEdit');
set(h,'String',num2str(sliderValue,'%3.4f'));


% --- Executes during object creation, after setting all properties.
function IsosurfaceSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IsosurfaceSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function IsosurfaceEdit_Callback(hObject, eventdata, handles)
% hObject    handle to IsosurfaceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IsosurfaceEdit as text
%        str2double(get(hObject,'String')) returns contents of IsosurfaceEdit as a double

editValue=str2num(get(handles.IsosurfaceEdit,'String'));
if (editValue>get(handles.IsosurfaceSlider,'Max'))
    editValue=get(handles.IsosurfaceSlider,'Max');
    set(handles.IsosurfaceEdit,'String',num2str(editValue,'%3.4f'));
elseif (editValue<get(handles.IsosurfaceSlider,'Min'))
    editValue=get(handles.IsosurfaceSlider,'Min');
    set(handles.IsosurfaceEdit,'String',num2str(editValue,'%3.4f'));
end
set(handles.IsosurfaceSlider,'Value',editValue);



        






% --- Executes during object creation, after setting all properties.
function IsosurfaceEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IsosurfaceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotButton.
function PlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filenameParameters=get(handles.SystemParametersFileEdit,'String');
filenameSysParams=setSysParamsFilename(filenameParameters);
load(filenameSysParams);
j=get(handles.ImaginaryTimesPopUpMenu,'Value');
filenamePsis=SysParams__filenames(j,:);

% Loading Psi values for plotting.
if(get(handles.LineRadioButton,'Value'))
    Npoints=200;
    [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,x,y,z]=...
        getPsisFull(filenamePsis,filenameSysParams,...
        get(handles.hedgehogizeCheckBox,'Value')*(1+get(handles.DiffFromHedgehogCheckebox,'Value')));
    theta=get(handles.LineThetaSlider,'Value');
    phi  =get(handles.LinePhiSlider,'Value');
    width=get(handles.LineWidthSlider,'Value');
    dlambda=width/(Npoints-1);
    x0=str2double(get(handles.LineX0Edit,'String'));
    y0=str2double(get(handles.LineY0Edit,'String'));
    z0=str2double(get(handles.LineZ0Edit,'String'));
    direction=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];
    r(1,:)=(-0.5*width+(0:(Npoints-1))*dlambda)*direction(1)+x0;
    r(2,:)=(-0.5*width+(0:(Npoints-1))*dlambda)*direction(2)+y0;
    r(3,:)=(-0.5*width+(0:(Npoints-1))*dlambda)*direction(3)+z0;
    xmin=-0.5*(SysParams__Mx-1)*SysParams__ax;
    ymin=-0.5*(SysParams__My-1)*SysParams__ay;
    zmin=-0.5*(SysParams__Mz-1)*SysParams__az;
    clear j;
    j(1,:)=floor((r(1,:)-xmin)/SysParams__ax)+1;
    j(2,:)=floor((r(2,:)-ymin)/SysParams__ay)+1;
    j(3,:)=floor((r(3,:)-zmin)/SysParams__az)+1;
    u(1,:)=r(1,:)-xmin-j(1,:)*SysParams__ax;
    u(2,:)=r(2,:)-ymin-j(2,:)*SysParams__ay;
    u(3,:)=r(3,:)-zmin-j(3,:)*SysParams__az;
    lineVector=linspace(-0.5*width,0.5*width,Npoints);
    for k=1:Npoints
        j(1,k)=min([max([j(1,k),1]),SysParams__Mx-1]);
        j(2,k)=min([max([j(2,k),1]),SysParams__My-1]);
        j(3,k)=min([max([j(3,k),1]),SysParams__Mz-1]);
    end
    tmpArray=Psi1;
    clear Psi1;
    for k=1:Npoints
        Psi1(k)=(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
               +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=Psi2;
    clear Psi2;
    for k=1:Npoints
        Psi2(k)=(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
               +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=Psi3;
    clear Psi3;
    for k=1:Npoints
        Psi3(k)=(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
               +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=n11;
    clear n11;
    for k=1:Npoints
        n11(k)  =(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
               +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=n12;
    clear n12;
    for k=1:Npoints
        n12(k)  =(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
               +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=n13;
    clear n13;
    for k=1:Npoints
        n13(k)  =(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
               +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=n22;
    clear n22;
    for k=1:Npoints
        n22(k)  =(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
               +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=n23;
    clear n23;
    for k=1:Npoints
        n23(k)  =(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
               +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=n33;
    clear n33;
    for k=1:Npoints
        n33(k)  =(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
               +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                 ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=nTOF11;
    clear nTOF11;
    for k=1:Npoints
        nTOF11(k)  =(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                    ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
                  +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                    ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=nTOF12;
    clear nTOF12;
    for k=1:Npoints
        nTOF12(k)  =(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                    ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
                  +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                    ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=nTOF13;
    clear nTOF13;
    for k=1:Npoints
        nTOF13(k)  =(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                    ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
                  +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                    ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=nTOF22;
    clear nTOF22;
    for k=1:Npoints
        nTOF22(k)  =(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                    ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
                  +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                    ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=nTOF23;
    clear nTOF23;
    for k=1:Npoints
        nTOF23(k)  =(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                    ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
                  +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                    ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
    tmpArray=nTOF33;
    clear nTOF33;
    for k=1:Npoints
        nTOF33(k)  =(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)  )).*(1-u(2,k))+...
                    ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)  )+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)  )).*u(2,k)    ).*(1-u(3,k))...
                  +(((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)  ,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)  ,j(3,k)+1)).*(1-u(2,k))+...
                    ((1-u(1,k)).*tmpArray(j(1,k)  ,j(2,k)+1,j(3,k)+1)+u(1,k).*tmpArray(j(1,k)+1,j(2,k)+1,j(3,k)+1)).*u(2,k)    ).*u(3,k);
    end
end
if(get(handles.FixedXRadioButton,'Value'))
    [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,y,z]=getPsisFixedX(filenamePsis,filenameSysParams, ...
          str2double(get(handles.FixedXEdit,'String')),...
          get(handles.hedgehogizeCheckBox,'Value')*(1+get(handles.DiffFromHedgehogCheckebox,'Value')));
end
if(get(handles.FixedYRadioButton,'Value'))
    [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,x,z]=getPsisFixedY(filenamePsis,filenameSysParams, ...
          str2double(get(handles.FixedYEdit,'String')),...
          get(handles.hedgehogizeCheckBox,'Value')*(1+get(handles.DiffFromHedgehogCheckebox,'Value')));
end
if(get(handles.FixedZRadioButton,'Value'))
    [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,x,y]=getPsisFixedZ(filenamePsis,filenameSysParams, ...
          str2double(get(handles.FixedZEdit,'String')),...
          get(handles.hedgehogizeCheckBox,'Value')*(1+get(handles.DiffFromHedgehogCheckebox,'Value')));
end
if(get(handles.IntegralRadioButton,'Value'))
    [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,x,y,z]=getPsisFull(filenamePsis,filenameSysParams,...
          get(handles.hedgehogizeCheckBox,'Value')*(1+get(handles.DiffFromHedgehogCheckebox,'Value')));
end
if(get(handles.EllipsoidRadioButton,'Value'))
    N_points=min([SysParams__Mx,SysParams__My,SysParams__Mz,200]);
    rx=str2num(get(handles.EllipsoidRxEdit,'String'));
    ry=str2num(get(handles.EllipsoidRyEdit,'String'));
    rz=str2num(get(handles.EllipsoidRzEdit,'String'));
    [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,X,Y,Z]=getPsisEllipsoid(filenamePsis,filenameSysParams,...
          rx,ry,rz,N_points,get(handles.hedgehogizeCheckBox,'Value')*(1+get(handles.DiffFromHedgehogCheckebox,'Value')));
end
if(get(handles.IsosurfaceRadioButton,'Value'))
    idealNoOfPoints=100;
    intersp_x=min([SysParams__Mx-1,max([1,round(SysParams__Mx/idealNoOfPoints)])]);
    intersp_y=min([SysParams__My-1,max([1,round(SysParams__My/idealNoOfPoints)])]);
    intersp_z=min([SysParams__Mz-1,max([1,round(SysParams__Mz/idealNoOfPoints)])]);
    [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,x,y,z]=getPsisInterspaced(filenamePsis,filenameSysParams,...
          intersp_x,intersp_y,intersp_z,get(handles.hedgehogizeCheckBox,'Value')*(1+get(handles.DiffFromHedgehogCheckebox,'Value')));
end

% Setting the array to plot
if get(handles.AbsPsiRadioButton,'Value')
    WhatToPlot=get(handles.AbsPsiRadioButton,'String');
elseif get(handles.ArgPsiRadioButton,'Value')
    WhatToPlot=get(handles.ArgPsiRadioButton,'String');
elseif get(handles.RePsiRadioButton,'Value')
    WhatToPlot=get(handles.RePsiRadioButton,'String');
elseif get(handles.ImPsiRadioButton,'Value')
    WhatToPlot=get(handles.ImPsiRadioButton,'String');
elseif get(handles.ArgPsiRadioButton,'Value')
    WhatToPlot=get(handles.ArgPsiRadioButton,'String');
elseif get(handles.nRadioButton,'Value')
    WhatToPlot=get(handles.nRadioButton,'String');
elseif get(handles.nTOFRadioButton,'Value')
    WhatToPlot=get(handles.nTOFRadioButton,'String');
elseif get(handles.AbsPsi_2RadioButton,'Value')
    WhatToPlot=get(handles.AbsPsi_2RadioButton,'String');
end
WhatToPlot=strtrim(WhatToPlot);
ComponentNo=get(handles.WhatToPlotPopUpMenu,'Value');
Component=get(handles.WhatToPlotPopUpMenu,'String');
Component=char(Component(ComponentNo));
if get(handles.HedgehogBasisRadioButton,'Value')==1
   whichBasis=1;
else
   whichBasis=2;
end
ArrayToPlot=setArrayToPlot(Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,WhatToPlot,Component,whichBasis);

% Plotting
if(get(handles.LineRadioButton,'Value'))
    if (strcmp(strtrim(Component),'1') || ...
        strcmp(strtrim(Component),'2') || ...
        strcmp(strtrim(Component),'3'))
        ArrayToPlot1=setArrayToPlot(Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,WhatToPlot,'1',whichBasis);
        ArrayToPlot2=setArrayToPlot(Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,WhatToPlot,'2',whichBasis);
        ArrayToPlot3=setArrayToPlot(Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,WhatToPlot,'3',whichBasis);
        figureTitle=[WhatToPlot,'1, ',WhatToPlot,'2, ',WhatToPlot,'3, ',...
                    ': Lineplot, theta=',num2str(theta),...
                    ', phi=',num2str(phi),', width=',num2str(width),...
                    ', x0=',num2str(x0),', y0=',num2str(y0),', z0=',num2str(z0)];
        plotLine3(lineVector,ArrayToPlot1,[WhatToPlot,'_1'],...
                   ArrayToPlot2,[WhatToPlot,'_2'],...
                   ArrayToPlot3,[WhatToPlot,'_3'],...
                   1,figureTitle)
    elseif (strcmp(strtrim(Component),'1-2') || ...
            strcmp(strtrim(Component),'2-3') || ...
            strcmp(strtrim(Component),'3-1'))
        ArrayToPlot1=setArrayToPlot(Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,WhatToPlot,'1-2',whichBasis);
        ArrayToPlot2=setArrayToPlot(Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,WhatToPlot,'2-3',whichBasis);
        ArrayToPlot3=setArrayToPlot(Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,WhatToPlot,'3-1',whichBasis);
        figureTitle=[WhatToPlot,'1-2, ',WhatToPlot,'2-3, ',WhatToPlot,'3-1, ',...
                    ': Lineplot, theta=',num2str(theta),...
                    ', phi=',num2str(phi),', width=',num2str(width),...
                    ', x0=',num2str(x0),', y0=',num2str(y0),', z0=',num2str(z0)];
        plotLine3(lineVector,ArrayToPlot1,[WhatToPlot,'1-2'],...
                   ArrayToPlot2,[WhatToPlot,'2-3'],...
                   ArrayToPlot3,[WhatToPlot,'3-1'],...
                   1,figureTitle)
    else
        figureTitle=[WhatToPlot,' ',Component,...
                    ': Lineplot, theta=',num2str(theta),...
                    ', phi=',num2str(phi),', width=',num2str(width),...
                    ', x0=',num2str(x0),', y0=',num2str(y0),', z0=',num2str(z0)];
        plotLine1(lineVector,ArrayToPlot,1,figureTitle);
    end
end
if(get(handles.FixedXRadioButton,'Value'))
    colorMap='default';
    figureNum=1;
    figureTitle=[WhatToPlot,', ',Component,', ',': Fixed x=',...
                 get(handles.FixedXEdit,'String')];
    if get(handles.FixedX2D3DRadioButton_2D,'Value')
        plotFlat(y,z,ArrayToPlot,2,figureNum,colorMap,figureTitle,'y','z');
    elseif get(handles.FixedX2D3DRadioButton_3D,'Value')
        plotFlat(y,z,ArrayToPlot,3,figureNum,colorMap,figureTitle,'y','z');
    end
end
if(get(handles.FixedYRadioButton,'Value'))
    colorMap='default';
    figureNum=1;
    figureTitle=[WhatToPlot,', ',Component,', ',': Fixed y=',...
                 get(handles.FixedYEdit,'String')];
    if get(handles.FixedY2D3DRadioButton_2D,'Value')
        plotFlat(x,z,ArrayToPlot,2,figureNum,colorMap,figureTitle,'x','z');
    elseif get(handles.FixedY2D3DRadioButton_3D,'Value')
        plotFlat(x,z,ArrayToPlot,3,figureNum,colorMap,figureTitle,'x','z');
    end
end
if(get(handles.FixedZRadioButton,'Value'))
    colorMap='default';
    figureNum=1;
    figureTitle=[WhatToPlot,', ',Component,', ',': Fixed z=',...
                 get(handles.FixedZEdit,'String')];
    if get(handles.FixedZ2D3DRadioButton_2D,'Value')
        plotFlat(x,y,ArrayToPlot,2,figureNum,colorMap,figureTitle,'x','y');
    elseif get(handles.FixedZ2D3DRadioButton_3D,'Value')
        plotFlat(x,y,ArrayToPlot,3,figureNum,colorMap,figureTitle,'x','y');
    end
end
if(get(handles.IntegralRadioButton,'Value'))
    colorMap='default';
    figureNum=1;
    if get(handles.Integral2D3DRadioButton_2D,'Value')
        plotDim=2;
    elseif get(handles.Integral2D3DRadioButton_3D,'Value')
        plotDim=3;
    end
    figureTitle=[WhatToPlot,' ',Component,': Integral, theta=',...
                 get(handles.IntegralThetaEdit,'String'),...
                 ', phi=',get(handles.IntegralPhiEdit,'String')];
    if get(handles.IntegralXAxisRadioButton,'Value')
        figureTitle=[figureTitle,' (x axis)'];
        plotIntegralFixedAxis(1,ArrayToPlot,plotDim,figureNum,colorMap,...
                              figureTitle,filenameSysParams);
    elseif get(handles.IntegralYAxisRadioButton,'Value')
        figureTitle=[figureTitle,' (y axis)'];
        plotIntegralFixedAxis(2,ArrayToPlot,plotDim,figureNum,colorMap,...
                              figureTitle,filenameSysParams);
    elseif get(handles.IntegralZAxisRadioButton,'Value')
        figureTitle=[figureTitle,' (z axis)'];
        plotIntegralFixedAxis(3,ArrayToPlot,plotDim,figureNum,colorMap,...
                              figureTitle,filenameSysParams);
    else
        theta=get(handles.IntegralThetaSlider,'Value');
        phi=get(handles.IntegralPhiSlider,'Value');
        plotIntegral(theta,phi,ArrayToPlot,plotDim,figureNum,colorMap,...
                     figureTitle,filenameSysParams)
    end
end
if(get(handles.EllipsoidRadioButton,'Value'))
    colorMap='default';
    figureNum=1;
    figureTitle=[WhatToPlot,', ',Component,', ',': Ellipsoid(rx=',...
                 get(handles.EllipsoidRxEdit,'String'),...
                 ', ry=',get(handles.EllipsoidRyEdit,'String'),...
                 ', rz=',get(handles.EllipsoidRzEdit,'String'),')'];
    plotEllipsoid(X,Y,Z,ArrayToPlot,figureNum,colorMap,...
                  figureTitle,'x','y','z');
    if(get(handles.EllipsoidLightCheckBox,'Value'))
        camlight;
        lighting 'phong';
    end
end
if(get(handles.IsosurfaceRadioButton,'Value'))
    isoValue=str2double(get(handles.IsosurfaceEdit,'String'));
    figureNum=1;
    figureTitle=['Isosurface: ',strtrim(WhatToPlot),...
                 ', ',strtrim(Component),'=',...
                 get(handles.IsosurfaceEdit,'String')];
    plotIsosurface(x,y,z,ArrayToPlot,isoValue, ...
                   figureNum,figureTitle,'x','y','z');
end
clear ArrayToPlot;




% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Create filename
filenameParameters=get(handles.SystemParametersFileEdit,'String');
filenameSysParams=setSysParamsFilename(filenameParameters);
load(filenameSysParams);
j=get(handles.ImaginaryTimesPopUpMenu,'Value');
filenameSave=regexprep(SysParams__filenames(j,:),' ','_');
filenameSave=[filenameSave,'___'];


if(get(handles.LineRadioButton,'Value'))
    filenameSave=[filenameSave,'Line(theta=',...
                  get(handles.LineThetaEdit,'String'),...
                  '_phi=',get(handles.LinePhiEdit,'String'),...
                  '_width=',get(handles.LineWidthEdit,'String'),...
                  '_x0=',get(handles.LineX0Edit,'String'),...
                  '_y0=',get(handles.LineY0Edit,'String'),...
                  '_z0=',get(handles.LineZ0Edit,'String'),...
                  ')__'];
end
if(get(handles.FixedXRadioButton,'Value'))
    filenameSave=[filenameSave,'FixedX(x=',get(handles.FixedXEdit,'String'),')__'];
end
if(get(handles.FixedYRadioButton,'Value'))
    filenameSave=[filenameSave,'FixedY(y=',get(handles.FixedYEdit,'String'),')__'];
end
if(get(handles.FixedZRadioButton,'Value'))
    filenameSave=[filenameSave,'FixedZ(z=',get(handles.FixedZEdit,'String'),')__'];
end
if(get(handles.IntegralRadioButton,'Value'))
    filenameSave=[filenameSave,'Integral(theta=',...
                  get(handles.IntegralThetaEdit,'String'),...
                  '_phi=',get(handles.IntegralPhiEdit,'String'),...
                  ')__'];
end
if(get(handles.EllipsoidRadioButton,'Value'))
    filenameSave=[filenameSave,'Ellipsoid(rx=',...
                  get(handles.EllipsoidRxEdit,'String'),...
                  '_ry=',get(handles.EllipsoidRyEdit,'String'),...
                  '_rz=',get(handles.EllipsoidRzEdit,'String'),')__'];
end
if(get(handles.IsosurfaceRadioButton,'Value'))
    filenameSave=[filenameSave,'Isosurface__'];
end

if get(handles.AbsPsiRadioButton,'Value')
    WhatToPlot=get(handles.AbsPsiRadioButton,'String');
elseif get(handles.ArgPsiRadioButton,'Value')
    WhatToPlot=get(handles.ArgPsiRadioButton,'String');
elseif get(handles.RePsiRadioButton,'Value')
    WhatToPlot=get(handles.RePsiRadioButton,'String');
elseif get(handles.ImPsiRadioButton,'Value')
    WhatToPlot=get(handles.ImPsiRadioButton,'String');
elseif get(handles.ArgPsiRadioButton,'Value')
    WhatToPlot=get(handles.ArgPsiRadioButton,'String');
elseif get(handles.nRadioButton,'Value')
    WhatToPlot=get(handles.nRadioButton,'String');
elseif get(handles.nTOFRadioButton,'Value')
    WhatToPlot=get(handles.nTOFRadioButton,'String');
elseif get(handles.AbsPsi_2RadioButton,'Value')
    WhatToPlot=get(handles.AbsPsi_2RadioButton,'String');
end
WhatToPlot=strtrim(WhatToPlot);
ComponentNo=get(handles.WhatToPlotPopUpMenu,'Value');
Component=get(handles.WhatToPlotPopUpMenu,'String');
Component=char(Component(ComponentNo));
filenameSave=[filenameSave,WhatToPlot,'_',Component];

if(get(handles.IsosurfaceRadioButton,'Value'))
    filenameSave=[filenameSave,'=',get(handles.IsosurfaceEdit,'String')];
end


if(~isempty(strtrim(get(handles.SavingCommentEdit,'String'))))
    commentString=strtrim(get(handles.SavingCommentEdit,'String'));
    for j=1:length(commentString)
        if(commentString(j)==' ')
            commentString(j)='_';
        end
    end
    filenameSave=[filenameSave,'__',commentString];
end

% Save image
figureNum=1;
figure(figureNum);
if(get(handles.JPGFileFormatRadioButton,'Value'))
    fileFormat='jpg';
end
if(get(handles.EPSFileFormatRadioButton,'Value'))
    fileFormat='eps';
end
if(get(handles.PDFFileFormatRadioButton,'Value'))
    fileFormat='pdf';
end
if(get(handles.FIGFileFormatRadioButton,'Value'))
    fileFormat='fig';
end
saveas(gcf,[filenameSave,'.',fileFormat],fileFormat);




function PlottingOptionsButtonGroup_SelectionChangeFcn(hObject, eventdata)
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 

switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'LineRadioButton'
      %execute this code when 'Fixed x' radio button is selected
      h=findobj('Tag','IsosurfaceMinText');
      set(h,'String','Min');
      h=findobj('Tag','IsosurfaceMaxText');
      set(h,'String','Max');

    case 'FixedXRadioButton'
      %execute this code when 'Fixed x' radio button is selected
      h=findobj('Tag','IsosurfaceMinText');
      set(h,'String','Min');
      h=findobj('Tag','IsosurfaceMaxText');
      set(h,'String','Max');

    case 'FixedYRadioButton'
      %execute this code when 'Fixed y' radio button is selected
      h=findobj('Tag','IsosurfaceMinText');
      set(h,'String','Min');
      h=findobj('Tag','IsosurfaceMaxText');
      set(h,'String','Max');

    case 'FixedZRadioButton'
      %execute this code when 'Fixed z' radio button is selected
      h=findobj('Tag','IsosurfaceMinText');
      set(h,'String','Min');
      h=findobj('Tag','IsosurfaceMaxText');
      set(h,'String','Max');
      
    case 'IntegralRadioButton'
      %execute this code when 'Fixed x' radio button is selected
      h=findobj('Tag','IsosurfaceMinText');
      set(h,'String','Min');
      h=findobj('Tag','IsosurfaceMaxText');
      set(h,'String','Max');

    case 'EllipsoidRadioButton'
      %execute this code when 'Ellipsoid' radio button is selected
      h=findobj('Tag','IsosurfaceMinText');
      set(h,'String','Min');
      h=findobj('Tag','IsosurfaceMaxText');
      set(h,'String','Max');
      
    case 'IsosurfaceRadioButton'
      %execute this code when 'Isosurface' radio button is selected
      %Setting Min and Max value and the sidebar
            h=findobj('Tag','SystemParametersFileEdit');
            filenameParameters=get(h,'String');
            filenameSysParams=setSysParamsFilename(filenameParameters);
            load(filenameSysParams);
            h=findobj('Tag','ImaginaryTimesPopUpMenu');
            j=get(h,'Value');
            filenamePsis=SysParams__filenames(j,:);
            % Get a low resolution sample of the data, to estimate the minimal 
            % and maximal values.
            intersp_x=min([SysParams__Mx,max([1,round(SysParams__Mx/64)])]);
            intersp_y=min([SysParams__My,max([1,round(SysParams__My/64)])]);
            intersp_z=min([SysParams__Mz,max([1,round(SysParams__Mz/64)])]);
            [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33]=getPsisInterspaced(filenamePsis,filenameSysParams,...
                  intersp_x,intersp_y,intersp_z,get(handles.hedgehogizeCheckBox,'Value')*(1+get(handles.DiffFromHedgehogCheckebox,'Value')));
            if get(handles.AbsPsiRadioButton,'Value')
                WhatToPlot=get(handles.AbsPsiRadioButton,'String');
            elseif get(handles.ArgPsiRadioButton,'Value')
                WhatToPlot=get(handles.ArgPsiRadioButton,'String');
            elseif get(handles.RePsiRadioButton,'Value')
                WhatToPlot=get(handles.RePsiRadioButton,'String');
            elseif get(handles.ImPsiRadioButton,'Value')
                WhatToPlot=get(handles.ImPsiRadioButton,'String');
            elseif get(handles.ArgPsiRadioButton,'Value')
                WhatToPlot=get(handles.ArgPsiRadioButton,'String');
            elseif get(handles.nRadioButton,'Value')
                WhatToPlot=get(handles.nRadioButton,'String');
            elseif get(handles.nTOFRadioButton,'Value')
                WhatToPlot=get(handles.nTOFRadioButton,'String');
            elseif get(handles.AbsPsi_2RadioButton,'Value')
                WhatToPlot=get(handles.AbsPsi_2RadioButton,'String');
            end
            WhatToPlot=strtrim(WhatToPlot);
            ComponentNo=get(handles.WhatToPlotPopUpMenu,'Value');
            Component=get(handles.WhatToPlotPopUpMenu,'String');
            Component=char(Component(ComponentNo));
            
            if get(handles.HedgehogBasisRadioButton,'Value')==1
                whichBasis=1;
            else
                whichBasis=2;
            end
            ArrayToPlot=setArrayToPlot(Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,WhatToPlot,Component,whichBasis);
            s=size(ArrayToPlot);
            jmax=size(s); jmax=jmax(2);
            minValue=min(ArrayToPlot);
            maxValue=max(ArrayToPlot);
            for j=2:jmax
                minValue=min(minValue);
                maxValue=max(maxValue);
            end
            if(minValue==maxValue)
                minValue=minValue-1e-5*abs(minValue);
                maxValue=maxValue+1e-5*abs(maxValue);
                if(minValue==0)
                    minValue=-1e-5;
                    maxValue=1e-5;
                end
            end
    
            clear  ArrayToPlot filenameParameters ...
                   filenameSysParams j s jmax dataTypeName ...
                   intersp_x intersp_y intersp_z;
            h=findobj('Tag','IsosurfaceMinText');
            set(h,'String',num2str(minValue,'%3.4f'));
            h=findobj('Tag','IsosurfaceMaxText');
            set(h,'String',num2str(maxValue,'%3.4f'));
            h=findobj('Tag','IsosurfaceSlider');
            set(h,'Min',minValue);
            set(h,'Max',maxValue);
            set(h,'Value',0.5*(minValue+maxValue));
            set(handles.IsosurfaceEdit,'String',...
                num2str(0.5*(minValue+maxValue),'%3.4f'));
    otherwise
       % Code for when there is no match.

end
%updates the handles structure
guidata(hObject, handles);


function WhatToPlotButtonGroup_SelectionChangeFcn(hObject, eventdata)

%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 

if(get(handles.IsosurfaceRadioButton,'Value'))
    % The slider and the texts have to be set when the
    % IsosurfaceRadioButton is on.
          h=findobj('Tag','SystemParametersFileEdit');
          filenameParameters=get(h,'String');
          filenameSysParams=setSysParamsFilename(filenameParameters);
          load(filenameSysParams);
          h=findobj('Tag','ImaginaryTimesPopUpMenu');
          j=get(h,'Value');
          filenamePsis=SysParams__filenames(j,:);
          % Get a low resolution sample of the data, to estimate the minimal 
          % and maximal values.
          intersp_x=min([SysParams__Mx,max([1,round(SysParams__Mx/64)])]);
          intersp_y=min([SysParams__My,max([1,round(SysParams__My/64)])]);
          intersp_z=min([SysParams__Mz,max([1,round(SysParams__Mz/64)])]);
          [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33]=getPsisInterspaced(filenamePsis,filenameSysParams,...
                intersp_x,intersp_y,intersp_z,get(handles.hedgehogizeCheckBox,'Value')*(1+get(handles.DiffFromHedgehogCheckebox,'Value')));
          % Chose which datatype is selected in WhatToPlotButtonGroup
          WhatToPlot=strtrim(get(eventdata.NewValue,'String'));
          ComponentNo=get(handles.WhatToPlotPopUpMenu,'Value');
          Component=get(handles.WhatToPlotPopUpMenu,'String');
          Component=char(Component(ComponentNo));
          if get(handles.HedgehogBasisRadioButton,'Value')==1
              whichBasis=1;
          else
              whichBasis=2;
          end
          ArrayToPlot=setArrayToPlot(Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,WhatToPlot,Component,whichBasis);
          s=size(ArrayToPlot);
          jmax=size(s); jmax=jmax(2);
          minValue=min(ArrayToPlot);
          maxValue=max(ArrayToPlot);
          for j=2:jmax
              minValue=min(minValue);
              maxValue=max(maxValue);
          end
          if(minValue==maxValue)
              minValue=minValue-1e-5*abs(minValue);
              maxValue=maxValue+1e-5*abs(maxValue);
              if(minValue==0)
                  minValue=-1e-5;
                  maxValue=1e-5;
              end
          end
          clear  ArrayToPlot filenameParameters filenameSysParams ...
                 j s jmax dataTypeName intersp_x intersp_y intersp_z;
          h=findobj('Tag','IsosurfaceMinText');
          set(h,'String',num2str(minValue,'%3.4f'));
          h=findobj('Tag','IsosurfaceMaxText');
          set(h,'String',num2str(maxValue,'%3.4f'));
          h=findobj('Tag','IsosurfaceSlider');
          set(h,'Min',minValue);
          set(h,'Max',maxValue);
          set(h,'Value',0.5*(minValue+maxValue));
          h=findobj('Tag','IsosurfaceEdit');
          set(h,'String',num2str(0.5*(minValue+maxValue),'%3.4f'));
end

%updates the handles structure
guidata(hObject, handles);


% --- Executes on button press in EllipsoidLightCheckBox.
function EllipsoidLightCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to EllipsoidLightCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EllipsoidLightCheckBox



function SavingCommentEdit_Callback(hObject, eventdata, handles)
% hObject    handle to SavingCommentEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SavingCommentEdit as text
%        str2double(get(hObject,'String')) returns contents of SavingCommentEdit as a double


% --- Executes during object creation, after setting all properties.
function SavingCommentEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SavingCommentEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in WhatToPlotPopUpMenu.
function WhatToPlotPopUpMenu_Callback(hObject, eventdata, handles)
% hObject    handle to WhatToPlotPopUpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns WhatToPlotPopUpMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from WhatToPlotPopUpMenu

if(get(handles.IsosurfaceRadioButton,'Value'))
    % The slider and the texts have to be set when the
    % IsosurfaceRadioButton is on.
          h=findobj('Tag','SystemParametersFileEdit');
          filenameParameters=get(h,'String');
          filenameSysParams=setSysParamsFilename(filenameParameters);
          load(filenameSysParams);
          h=findobj('Tag','ImaginaryTimesPopUpMenu');
          j=get(h,'Value');
          filenamePsis=SysParams__filenames(j,:);
          % Get a low resolution sample of the data, to estimate the minimal 
          % and maximal values.
          intersp_x=min([SysParams__Mx,max([1,round(SysParams__Mx/64)])]);
          intersp_y=min([SysParams__My,max([1,round(SysParams__My/64)])]);
          intersp_z=min([SysParams__Mz,max([1,round(SysParams__Mz/64)])]);
          [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33]=getPsisInterspaced(filenamePsis,filenameSysParams,...
                intersp_x,intersp_y,intersp_z,get(handles.hedgehogizeCheckBox,'Value')*(1+get(handles.DiffFromHedgehogCheckebox,'Value')));
          % Chose which datatype is selected in WhatToPlotButtonGroup
          if(get(handles.AbsPsiRadioButton,'Value'))
              WhatToPlot=get(handles.AbsPsiRadioButton,'String');
          elseif(get(handles.ArgPsiRadioButton,'Value'))
              WhatToPlot=get(handles.ArgPsiRadioButton,'String');
          elseif(get(handles.RePsiRadioButton,'Value'))
              WhatToPlot=get(handles.RePsiRadioButton,'String');
          elseif(get(handles.ImPsiRadioButton,'Value'))
              WhatToPlot=get(handles.ImPsiRadioButton,'String');
          elseif(get(handles.nRadioButton,'Value'))
              WhatToPlot=get(handles.nRadioButton,'String');
          elseif(get(handles.nTOFRadioButton,'Value'))
              WhatToPlot=get(handles.nTOFRadioButton,'String');
          elseif(get(handles.AbsPsi_2RadioButton,'Value'))
              WhatToPlot=get(handles.AbsPsi_2RadioButton,'String');
          end
          WhatToPlot=strtrim(WhatToPlot);
          ComponentNo=get(handles.WhatToPlotPopUpMenu,'Value');
          Component=get(handles.WhatToPlotPopUpMenu,'String');
          Component=char(Component(ComponentNo));
          if get(handles.HedgehogBasisRadioButton,'Value')==1
              whichBasis=1;
          else
              whichBasis=2;
          end
          ArrayToPlot=setArrayToPlot(Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,...
                                    nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,...
                                    WhatToPlot,Component,whichBasis);
          s=size(ArrayToPlot);
          jmax=size(s); jmax=jmax(2);
          minValue=min(ArrayToPlot);
          maxValue=max(ArrayToPlot);
          for j=2:jmax
              minValue=min(minValue);
              maxValue=max(maxValue);
          end
          if(minValue==maxValue)
              minValue=minValue-1e-5*abs(minValue);
              maxValue=maxValue+1e-5*abs(maxValue);
              if(minValue==0)
                  minValue=-1e-5;
                  maxValue=1e-5;
              end
          end
          clear ArrayToPlot filenameParameters filenameSysParams ...
                j s jmax dataTypeName intersp_x intersp_y intersp_z;
          h=findobj('Tag','IsosurfaceMinText');
          set(h,'String',num2str(minValue,'%3.4f'));
          h=findobj('Tag','IsosurfaceMaxText');
          set(h,'String',num2str(maxValue,'%3.4f'));
          h=findobj('Tag','IsosurfaceSlider');
          set(h,'Min',minValue);
          set(h,'Max',maxValue);
          set(h,'Value',0.5*(minValue+maxValue));
          h=findobj('Tag','IsosurfaceEdit');
          set(h,'String',num2str(0.5*(minValue+maxValue),'%3.4f'));
end

%updates the handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function WhatToPlotPopUpMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhatToPlotPopUpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function IntegralThetaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to IntegralThetaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if get(handles.IntegralXAxisRadioButton,'Value')
    set(handles.IntegralThetaSlider,'Value',90.0);
    set(handles.IntegralThetaEdit,'String',num2str(90.0));
    set(handles.IntegralPhiSlider,'Value',0.0);
    set(handles.IntegralPhiEdit,'String',num2str(0.0));
elseif get(handles.IntegralYAxisRadioButton,'Value')
    set(handles.IntegralThetaSlider,'Value',90.0);
    set(handles.IntegralThetaEdit,'String',num2str(90.0));
    set(handles.IntegralPhiSlider,'Value',90.0);
    set(handles.IntegralPhiEdit,'String',num2str(90.0));
elseif get(handles.IntegralZAxisRadioButton,'Value')
    set(handles.IntegralThetaSlider,'Value',0.0);
    set(handles.IntegralThetaEdit,'String',num2str(0.0));
    set(handles.IntegralPhiSlider,'Value',0.0);
    set(handles.IntegralPhiEdit,'String',num2str(0.0));
else
    sliderValue=get(handles.IntegralThetaSlider,'Value');
    set(handles.IntegralThetaEdit,'String',num2str(sliderValue,'%3.4f'));
end



% --- Executes during object creation, after setting all properties.
function IntegralThetaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IntegralThetaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function IntegralThetaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to IntegralThetaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IntegralThetaEdit as text
%        str2double(get(hObject,'String')) returns contents of IntegralThetaEdit as a double

if get(handles.IntegralXAxisRadioButton,'Value')
    set(handles.IntegralThetaSlider,'Value',90.0);
    set(handles.IntegralThetaEdit,'String',num2str(90.0));
    set(handles.IntegralPhiSlider,'Value',0.0);
    set(handles.IntegralPhiEdit,'String',num2str(0.0));
elseif get(handles.IntegralYAxisRadioButton,'Value')
    set(handles.IntegralThetaSlider,'Value',90.0);
    set(handles.IntegralThetaEdit,'String',num2str(90.0));
    set(handles.IntegralPhiSlider,'Value',90.0);
    set(handles.IntegralPhiEdit,'String',num2str(90.0));
elseif get(handles.IntegralZAxisRadioButton,'Value')
    set(handles.IntegralThetaSlider,'Value',0.0);
    set(handles.IntegralThetaEdit,'String',num2str(0.0));
    set(handles.IntegralPhiSlider,'Value',0.0);
    set(handles.IntegralPhiEdit,'String',num2str(0.0));
else
    editValue=str2double(get(handles.IntegralThetaEdit,'String'));
    if (editValue>get(handles.IntegralThetaSlider,'Max'))
        editValue=get(handles.IntegralThetaSlider,'Max');
        set(handles.IntegralThetaEdit,'String',num2str(editValue,'%3.4f'));
    elseif (editValue<get(handles.IntegralThetaSlider,'Min'))
        editValue=get(handles.IntegralThetaSlider,'Min');
        set(handles.IntegralThetaEdit,'String',num2str(editValue,'%3.4f'));
    end
    set(handles.IntegralThetaSlider,'Value',editValue);
end



% --- Executes during object creation, after setting all properties.
function IntegralThetaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IntegralThetaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function IntegralPhiSlider_Callback(hObject, eventdata, handles)
% hObject    handle to IntegralPhiSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if get(handles.IntegralXAxisRadioButton,'Value')
    set(handles.IntegralThetaSlider,'Value',90.0);
    set(handles.IntegralThetaEdit,'String',num2str(90.0));
    set(handles.IntegralPhiSlider,'Value',0.0);
    set(handles.IntegralPhiEdit,'String',num2str(0.0));
elseif get(handles.IntegralYAxisRadioButton,'Value')
    set(handles.IntegralThetaSlider,'Value',90.0);
    set(handles.IntegralThetaEdit,'String',num2str(90.0));
    set(handles.IntegralPhiSlider,'Value',90.0);
    set(handles.IntegralPhiEdit,'String',num2str(90.0));
elseif get(handles.IntegralZAxisRadioButton,'Value')
    set(handles.IntegralThetaSlider,'Value',0.0);
    set(handles.IntegralThetaEdit,'String',num2str(0.0));
    set(handles.IntegralPhiSlider,'Value',0.0);
    set(handles.IntegralPhiEdit,'String',num2str(0.0));
else
    sliderValue=get(handles.IntegralPhiSlider,'Value');
    set(handles.IntegralPhiEdit,'String',num2str(sliderValue,'%3.4f'));
end



% --- Executes during object creation, after setting all properties.
function IntegralPhiSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IntegralPhiSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function IntegralPhiEdit_Callback(hObject, eventdata, handles)
% hObject    handle to IntegralPhiEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IntegralPhiEdit as text
%        str2double(get(hObject,'String')) returns contents of IntegralPhiEdit as a double

if get(handles.IntegralXAxisRadioButton,'Value')
    set(handles.IntegralThetaSlider,'Value',90.0);
    set(handles.IntegralThetaEdit,'String',num2str(90.0));
    set(handles.IntegralPhiSlider,'Value',0.0);
    set(handles.IntegralPhiEdit,'String',num2str(0.0));
elseif get(handles.IntegralYAxisRadioButton,'Value')
    set(handles.IntegralThetaSlider,'Value',90.0);
    set(handles.IntegralThetaEdit,'String',num2str(90.0));
    set(handles.IntegralPhiSlider,'Value',90.0);
    set(handles.IntegralPhiEdit,'String',num2str(90.0));
elseif get(handles.IntegralZAxisRadioButton,'Value')
    set(handles.IntegralThetaSlider,'Value',0.0);
    set(handles.IntegralThetaEdit,'String',num2str(0.0));
    set(handles.IntegralPhiSlider,'Value',0.0);
    set(handles.IntegralPhiEdit,'String',num2str(0.0));
else
    editValue=str2double(get(handles.IntegralPhiEdit,'String'));
    if (editValue>get(handles.IntegralPhiSlider,'Max'))
        editValue=get(handles.IntegralPhiSlider,'Max');
        set(handles.IntegralPhiEdit,'String',num2str(editValue,'%3.4f'));
    elseif (editValue<get(handles.IntegralPhiSlider,'Min'))
        editValue=get(handles.IntegralPhiSlider,'Min');
        set(handles.IntegralPhiEdit,'String',num2str(editValue,'%3.4f'));
    end
    set(handles.IntegralPhiSlider,'Value',editValue);
end



% --- Executes during object creation, after setting all properties.
function IntegralPhiEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IntegralPhiEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function LineThetaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to LineThetaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderValue=get(handles.LineThetaSlider,'Value');
set(handles.LineThetaEdit,'String',num2str(sliderValue,'%3.4f'));

filenameParameters=get(handles.SystemParametersFileEdit,'String');
filenameSysParams=setSysParamsFilename(filenameParameters);
load(filenameSysParams);
theta=get(handles.LineThetaSlider,'Value')/180*pi;
phi=get(handles.LinePhiSlider,'Value')/180*pi;
x0=str2double(get(handles.LineX0Edit,'String'));
y0=str2double(get(handles.LineY0Edit,'String'));
z0=str2double(get(handles.LineZ0Edit,'String'));
xmax= 0.5*(SysParams__Mx-1)*SysParams__ax;
xmin=-0.5*(SysParams__Mx-1)*SysParams__ax;
ymax= 0.5*(SysParams__My-1)*SysParams__ay;
ymin=-0.5*(SysParams__My-1)*SysParams__ay;
zmax= 0.5*(SysParams__Mz-1)*SysParams__az;
zmin=-0.5*(SysParams__Mz-1)*SysParams__az;

widthMax=sqrt((xmax-xmin)^2+(ymax-ymin)^2+(zmax-zmin)^2);
if (sin(theta)*cos(phi))~=0.0
    widthMax=min([abs((xmax-x0)/(sin(theta)*cos(phi))),...
                  abs((xmin-x0)/(sin(theta)*cos(phi)))]);
end
if (sin(theta)*sin(phi))~=0.0
    widthMax=min([widthMax,...
                  abs((ymax-y0)/(sin(theta)*sin(phi))),...
                  abs((ymin-y0)/(sin(theta)*sin(phi)))]);
end
if cos(theta)~=0.0
    widthMax=min([widthMax,...
                  abs((zmax-z0)/cos(theta)),...
                  abs((zmin-z0)/cos(theta))]);
end
widthMax=widthMax*2.0;
width=str2double(get(handles.LineWidthEdit,'String'));
width=min([width,widthMax]);
set(handles.LineWidthSlider,'Max',widthMax);
set(handles.LineWidthSlider,'Value',width);
set(handles.LineWidthMaxWidthText,'String',num2str(widthMax,'%3.3f'));
set(handles.LineWidthEdit,'String',num2str(width,'%3.3f'));




% --- Executes during object creation, after setting all properties.
function LineThetaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LineThetaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function LineThetaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to LineThetaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LineThetaEdit as text
%        str2double(get(hObject,'String')) returns contents of LineThetaEdit as a double

editValue=str2num(get(handles.LineThetaEdit,'String'));
if (editValue>get(handles.LineThetaSlider,'Max'))
    editValue=get(handles.LineThetaSlider,'Max');
    set(handles.LineThetaEdit,'String',num2str(editValue,'%3.4f'));
elseif (editValue<get(handles.LineThetaSlider,'Min'))
    editValue=get(handles.LineThetaSlider,'Min');
    set(handles.LineThetaEdit,'String',num2str(editValue,'%3.4f'));
end
set(handles.LineThetaSlider,'Value',editValue);


filenameParameters=get(handles.SystemParametersFileEdit,'String');
filenameSysParams=setSysParamsFilename(filenameParameters);
load(filenameSysParams);
theta=get(handles.LineThetaSlider,'Value')/180*pi;
phi=get(handles.LinePhiSlider,'Value')/180*pi;
x0=str2double(get(handles.LineX0Edit,'String'));
y0=str2double(get(handles.LineY0Edit,'String'));
z0=str2double(get(handles.LineZ0Edit,'String'));
xmax= 0.5*(SysParams__Mx-1)*SysParams__ax;
xmin=-0.5*(SysParams__Mx-1)*SysParams__ax;
ymax= 0.5*(SysParams__My-1)*SysParams__ay;
ymin=-0.5*(SysParams__My-1)*SysParams__ay;
zmax= 0.5*(SysParams__Mz-1)*SysParams__az;
zmin=-0.5*(SysParams__Mz-1)*SysParams__az;

widthMax=sqrt((xmax-xmin)^2+(ymax-ymin)^2+(zmax-zmin)^2);
if (sin(theta)*cos(phi))~=0.0
    widthMax=min([abs((xmax-x0)/(sin(theta)*cos(phi))),...
                  abs((xmin-x0)/(sin(theta)*cos(phi)))]);
end
if (sin(theta)*sin(phi))~=0.0
    widthMax=min([widthMax,...
                  abs((ymax-y0)/(sin(theta)*sin(phi))),...
                  abs((ymin-y0)/(sin(theta)*sin(phi)))]);
end
if cos(theta)~=0.0
    widthMax=min([widthMax,...
                  abs((zmax-z0)/cos(theta)),...
                  abs((zmin-z0)/cos(theta))]);
end
widthMax=widthMax*2.0;
width=str2double(get(handles.LineWidthEdit,'String'));
width=min([width,widthMax]);
set(handles.LineWidthSlider,'Max',widthMax);
set(handles.LineWidthSlider,'Value',width);
set(handles.LineWidthMaxWidthText,'String',num2str(widthMax,'%3.3f'));
set(handles.LineWidthEdit,'String',num2str(width,'%3.3f'));



% --- Executes during object creation, after setting all properties.
function LineThetaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LineThetaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function LinePhiSlider_Callback(hObject, eventdata, handles)
% hObject    handle to LinePhiSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderValue=get(handles.LinePhiSlider,'Value');
set(handles.LinePhiEdit,'String',num2str(sliderValue,'%3.4f'));


filenameParameters=get(handles.SystemParametersFileEdit,'String');
filenameSysParams=setSysParamsFilename(filenameParameters);
load(filenameSysParams);
theta=get(handles.LineThetaSlider,'Value')/180*pi;
phi=get(handles.LinePhiSlider,'Value')/180*pi;
x0=str2double(get(handles.LineX0Edit,'String'));
y0=str2double(get(handles.LineY0Edit,'String'));
z0=str2double(get(handles.LineZ0Edit,'String'));
xmax= 0.5*(SysParams__Mx-1)*SysParams__ax;
xmin=-0.5*(SysParams__Mx-1)*SysParams__ax;
ymax= 0.5*(SysParams__My-1)*SysParams__ay;
ymin=-0.5*(SysParams__My-1)*SysParams__ay;
zmax= 0.5*(SysParams__Mz-1)*SysParams__az;
zmin=-0.5*(SysParams__Mz-1)*SysParams__az;

widthMax=sqrt((xmax-xmin)^2+(ymax-ymin)^2+(zmax-zmin)^2)*0.5;
if (sin(theta)*cos(phi))~=0.0
    widthMax=min([abs((xmax-x0)/(sin(theta)*cos(phi))),...
                  abs((xmin-x0)/(sin(theta)*cos(phi)))]);
end
if (sin(theta)*sin(phi))~=0.0
    widthMax=min([widthMax,...
                  abs((ymax-y0)/(sin(theta)*sin(phi))),...
                  abs((ymin-y0)/(sin(theta)*sin(phi)))]);
end
if cos(theta)~=0.0
    widthMax=min([widthMax,...
                  abs((zmax-z0)/cos(theta)),...
                  abs((zmin-z0)/cos(theta))]);
end
widthMax=widthMax*2.0;
width=str2double(get(handles.LineWidthEdit,'String'));
width=min([width,widthMax]);
set(handles.LineWidthSlider,'Max',widthMax);
set(handles.LineWidthSlider,'Value',width);
set(handles.LineWidthMaxWidthText,'String',num2str(widthMax,'%3.3f'));
set(handles.LineWidthEdit,'String',num2str(width,'%3.3f'));



% --- Executes during object creation, after setting all properties.
function LinePhiSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LinePhiSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function LinePhiEdit_Callback(hObject, eventdata, handles)
% hObject    handle to LinePhiEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LinePhiEdit as text
%        str2double(get(hObject,'String')) returns contents of LinePhiEdit as a double

editValue=str2num(get(handles.LinePhiEdit,'String'));
if (editValue>get(handles.LinePhiSlider,'Max'))
    editValue=get(handles.LinePhiSlider,'Max');
    set(handles.LinePhiEdit,'String',num2str(editValue,'%3.4f'));
elseif (editValue<get(handles.LinePhiSlider,'Min'))
    editValue=get(handles.LinePhiSlider,'Min');
    set(handles.LinePhiEdit,'String',num2str(editValue,'%3.4f'));
end
set(handles.LinePhiSlider,'Value',editValue);


filenameParameters=get(handles.SystemParametersFileEdit,'String');
filenameSysParams=setSysParamsFilename(filenameParameters);
load(filenameSysParams);
theta=get(handles.LineThetaSlider,'Value')/180*pi;
phi=get(handles.LinePhiSlider,'Value')/180*pi;
x0=str2double(get(handles.LineX0Edit,'String'));
y0=str2double(get(handles.LineY0Edit,'String'));
z0=str2double(get(handles.LineZ0Edit,'String'));
xmax= 0.5*(SysParams__Mx-1)*SysParams__ax;
xmin=-0.5*(SysParams__Mx-1)*SysParams__ax;
ymax= 0.5*(SysParams__My-1)*SysParams__ay;
ymin=-0.5*(SysParams__My-1)*SysParams__ay;
zmax= 0.5*(SysParams__Mz-1)*SysParams__az;
zmin=-0.5*(SysParams__Mz-1)*SysParams__az;

widthMax=sqrt((xmax-xmin)^2+(ymax-ymin)^2+(zmax-zmin)^2);
if (sin(theta)*cos(phi))~=0.0
    widthMax=min([abs((xmax-x0)/(sin(theta)*cos(phi))),...
                  abs((xmin-x0)/(sin(theta)*cos(phi)))]);
end
if (sin(theta)*sin(phi))~=0.0
    widthMax=min([widthMax,...
                  abs((ymax-y0)/(sin(theta)*sin(phi))),...
                  abs((ymin-y0)/(sin(theta)*sin(phi)))]);
end
if cos(theta)~=0.0
    widthMax=min([widthMax,...
                  abs((zmax-z0)/cos(theta)),...
                  abs((zmin-z0)/cos(theta))]);
end
widthMax=widthMax*2.0;
width=str2double(get(handles.LineWidthEdit,'String'));
width=min([width,widthMax]);
set(handles.LineWidthSlider,'Max',widthMax);
set(handles.LineWidthSlider,'Value',width);
set(handles.LineWidthMaxWidthText,'String',num2str(widthMax,'%3.3f'));
set(handles.LineWidthEdit,'String',num2str(width,'%3.3f'));



% --- Executes during object creation, after setting all properties.
function LinePhiEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LinePhiEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function LineWidthSlider_Callback(hObject, eventdata, handles)
% hObject    handle to LineWidthSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderValue=get(handles.LineWidthSlider,'Value');
set(handles.LineWidthEdit,'String',num2str(sliderValue,'%3.4f'));



% --- Executes during object creation, after setting all properties.
function LineWidthSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LineWidthSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function LineWidthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to LineWidthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LineWidthEdit as text
%        str2double(get(hObject,'String')) returns contents of LineWidthEdit as a double

editValue=str2num(get(handles.LineWidthEdit,'String'));
if (editValue>get(handles.LineWidthSlider,'Max'))
    editValue=get(handles.LineWidthSlider,'Max');
    set(handles.LineWidthEdit,'String',num2str(editValue,'%3.4f'));
elseif (editValue<get(handles.LineWidthSlider,'Min'))
    editValue=get(handles.LineWidthSlider,'Min');
    set(handles.LineWidthEdit,'String',num2str(editValue,'%3.4f'));
end
set(handles.LineWidthSlider,'Value',editValue);


% --- Executes during object creation, after setting all properties.
function LineWidthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LineWidthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LineZ0Edit_Callback(hObject, eventdata, handles)
% hObject    handle to LineZ0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LineZ0Edit as text
%        str2double(get(hObject,'String')) returns contents of LineZ0Edit as a double
filenameParameters=get(handles.SystemParametersFileEdit,'String');
filenameSysParams=setSysParamsFilename(filenameParameters);
load(filenameSysParams);
theta=get(handles.LineThetaSlider,'Value')/180*pi;
phi=get(handles.LinePhiSlider,'Value')/180*pi;
z0=str2double(get(handles.LineZ0Edit,'String'));
z0=min([z0, 0.5*(SysParams__Mz-1)*SysParams__az*0.99]);
z0=max([z0,-0.5*(SysParams__Mz-1)*SysParams__az*0.99]);
set(handles.LineZ0Edit,'String',num2str(z0));
y0=str2double(get(handles.LineY0Edit,'String'));
x0=str2double(get(handles.LineX0Edit,'String'));
xmax= 0.5*(SysParams__Mx-1)*SysParams__ax;
xmin=-0.5*(SysParams__Mx-1)*SysParams__ax;
ymax= 0.5*(SysParams__My-1)*SysParams__ay;
ymin=-0.5*(SysParams__My-1)*SysParams__ay;
zmax= 0.5*(SysParams__Mz-1)*SysParams__az;
zmin=-0.5*(SysParams__Mz-1)*SysParams__az;

widthMax=sqrt((xmax-xmin)^2+(ymax-ymin)^2+(zmax-zmin)^2);
if (sin(theta)*cos(phi))~=0.0
    widthMax=min([abs((xmax-x0)/(sin(theta)*cos(phi))),...
                  abs((xmin-x0)/(sin(theta)*cos(phi)))]);
end
if (sin(theta)*sin(phi))~=0.0
    widthMax=min([widthMax,...
                  abs((ymax-y0)/(sin(theta)*sin(phi))),...
                  abs((ymin-y0)/(sin(theta)*sin(phi)))]);
end
if cos(theta)~=0.0
    widthMax=min([widthMax,...
                  abs((zmax-z0)/cos(theta)),...
                  abs((zmin-z0)/cos(theta))]);
end
widthMax=widthMax*2.0;
width=str2double(get(handles.LineWidthEdit,'String'));
width=min([width,widthMax]);
set(handles.LineWidthSlider,'Max',widthMax);
set(handles.LineWidthSlider,'Value',width);
set(handles.LineWidthMaxWidthText,'String',num2str(widthMax,'%3.3f'));
set(handles.LineWidthEdit,'String',num2str(width,'%3.3f'));



% --- Executes during object creation, after setting all properties.
function LineZ0Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LineZ0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LineY0Edit_Callback(hObject, eventdata, handles)
% hObject    handle to LineY0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LineY0Edit as text
%        str2double(get(hObject,'String')) returns contents of LineY0Edit as a double
filenameParameters=get(handles.SystemParametersFileEdit,'String');
filenameSysParams=setSysParamsFilename(filenameParameters);
load(filenameSysParams);
theta=get(handles.LineThetaSlider,'Value')/180*pi;
phi=get(handles.LinePhiSlider,'Value')/180*pi;
y0=str2double(get(handles.LineY0Edit,'String'));
y0=min([y0, 0.5*(SysParams__My-1)*SysParams__ay*0.99]);
y0=max([y0,-0.5*(SysParams__My-1)*SysParams__ay*0.99]);
set(handles.LineY0Edit,'String',num2str(y0));
x0=str2double(get(handles.LineX0Edit,'String'));
z0=str2double(get(handles.LineZ0Edit,'String'));
xmax= 0.5*(SysParams__Mx-1)*SysParams__ax;
xmin=-0.5*(SysParams__Mx-1)*SysParams__ax;
ymax= 0.5*(SysParams__My-1)*SysParams__ay;
ymin=-0.5*(SysParams__My-1)*SysParams__ay;
zmax= 0.5*(SysParams__Mz-1)*SysParams__az;
zmin=-0.5*(SysParams__Mz-1)*SysParams__az;

widthMax=sqrt((xmax-xmin)^2+(ymax-ymin)^2+(zmax-zmin)^2);
if (sin(theta)*cos(phi))~=0.0
    widthMax=min([abs((xmax-x0)/(sin(theta)*cos(phi))),...
                  abs((xmin-x0)/(sin(theta)*cos(phi)))]);
end
if (sin(theta)*sin(phi))~=0.0
    widthMax=min([widthMax,...
                  abs((ymax-y0)/(sin(theta)*sin(phi))),...
                  abs((ymin-y0)/(sin(theta)*sin(phi)))]);
end
if cos(theta)~=0.0
    widthMax=min([widthMax,...
                  abs((zmax-z0)/cos(theta)),...
                  abs((zmin-z0)/cos(theta))]);
end
widthMax=widthMax*2.0;
width=str2double(get(handles.LineWidthEdit,'String'));
width=min([width,widthMax]);
set(handles.LineWidthSlider,'Max',widthMax);
set(handles.LineWidthSlider,'Value',width);
set(handles.LineWidthMaxWidthText,'String',num2str(widthMax,'%3.3f'));
set(handles.LineWidthEdit,'String',num2str(width,'%3.3f'));



% --- Executes during object creation, after setting all properties.
function LineY0Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LineY0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LineX0Edit_Callback(hObject, eventdata, handles)
% hObject    handle to LineX0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LineX0Edit as text
%        str2double(get(hObject,'String')) returns contents of LineX0Edit as a double
filenameParameters=get(handles.SystemParametersFileEdit,'String');
filenameSysParams=setSysParamsFilename(filenameParameters);
load(filenameSysParams);
theta=get(handles.LineThetaSlider,'Value')/180*pi;
phi=get(handles.LinePhiSlider,'Value')/180*pi;
x0=str2double(get(handles.LineX0Edit,'String'));
x0=min([x0, 0.5*(SysParams__Mx-1)*SysParams__ax*0.99]);
x0=max([x0,-0.5*(SysParams__Mx-1)*SysParams__ax*0.99]);
set(handles.LineX0Edit,'String',num2str(x0));
y0=str2double(get(handles.LineY0Edit,'String'));
z0=str2double(get(handles.LineZ0Edit,'String'));
xmax= 0.5*(SysParams__Mx-1)*SysParams__ax;
xmin=-0.5*(SysParams__Mx-1)*SysParams__ax;
ymax= 0.5*(SysParams__My-1)*SysParams__ay;
ymin=-0.5*(SysParams__My-1)*SysParams__ay;
zmax= 0.5*(SysParams__Mz-1)*SysParams__az;
zmin=-0.5*(SysParams__Mz-1)*SysParams__az;

widthMax=sqrt((xmax-xmin)^2+(ymax-ymin)^2+(zmax-zmin)^2);
if (sin(theta)*cos(phi))~=0.0
    widthMax=min([abs((xmax-x0)/(sin(theta)*cos(phi))),...
                  abs((xmin-x0)/(sin(theta)*cos(phi)))]);
end
if (sin(theta)*sin(phi))~=0.0
    widthMax=min([widthMax,...
                  abs((ymax-y0)/(sin(theta)*sin(phi))),...
                  abs((ymin-y0)/(sin(theta)*sin(phi)))]);
end
if cos(theta)~=0.0
    widthMax=min([widthMax,...
                  abs((zmax-z0)/cos(theta)),...
                  abs((zmin-z0)/cos(theta))]);
end
widthMax=widthMax*2.0;
width=str2double(get(handles.LineWidthEdit,'String'));
width=min([width,widthMax]);
set(handles.LineWidthSlider,'Max',widthMax);
set(handles.LineWidthSlider,'Value',width);
set(handles.LineWidthMaxWidthText,'String',num2str(widthMax,'%3.3f'));
set(handles.LineWidthEdit,'String',num2str(width,'%3.3f'));



% --- Executes during object creation, after setting all properties.
function LineX0Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LineX0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FixedRatioCheckBox.
function FixedRatioCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to FixedRatioCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FixedRatioCheckBox

function IntegralFixedAxisButtonGroup_SelectionChangeFcn(hObject, eventdata)
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'IntegralXAxisRadioButton'
        h=findobj('Tag','IntegralThetaSlider'); set(h,'Value',90.0);
        h=findobj('Tag','IntegralThetaEdit');   set(h,'String',num2str(90.0));
        h=findobj('Tag','IntegralPhiSlider');   set(h,'Value',0.0);
        h=findobj('Tag','IntegralPhiEdit');     set(h,'String',num2str(0.0));
    case 'IntegralYAxisRadioButton'
        h=findobj('Tag','IntegralThetaSlider'); set(h,'Value',90.0);
        h=findobj('Tag','IntegralThetaEdit');   set(h,'String',num2str(90.0));
        h=findobj('Tag','IntegralPhiSlider');   set(h,'Value',90.0);
        h=findobj('Tag','IntegralPhiEdit');     set(h,'String',num2str(90.0));
    case 'IntegralZAxisRadioButton'
        h=findobj('Tag','IntegralThetaSlider'); set(h,'Value',0.0);
        h=findobj('Tag','IntegralThetaEdit');   set(h,'String',num2str(0.0));
        h=findobj('Tag','IntegralPhiSlider');   set(h,'Value',0.0);
        h=findobj('Tag','IntegralPhiEdit');     set(h,'String',num2str(0.0));
    otherwise
end


% --- Executes during object creation, after setting all properties.
function IntegralXAxisRadioButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IntegralXAxisRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in hedgehogizeCheckBox.
function hedgehogizeCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to hedgehogizeCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hedgehogizeCheckBox


% --- Executes on button press in DiffFromHedgehogCheckebox.
function DiffFromHedgehogCheckebox_Callback(hObject, eventdata, handles)
% hObject    handle to DiffFromHedgehogCheckebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DiffFromHedgehogCheckebox
