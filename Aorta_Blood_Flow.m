function varargout = Aorta_Blood_Flow(varargin)
% AORTA_BLOOD_FLOW MATLAB code for Aorta_Blood_Flow.fig
%      AORTA_BLOOD_FLOW, by itself, creates a new AORTA_BLOOD_FLOW or raises the existing
%      singleton*.
%
%      H = AORTA_BLOOD_FLOW returns the handle to a new AORTA_BLOOD_FLOW or the handle to
%      the existing singleton*.
%
%      AORTA_BLOOD_FLOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AORTA_BLOOD_FLOW.M with the given input arguments.
%
%      AORTA_BLOOD_FLOW('Property','Value',...) creates a new AORTA_BLOOD_FLOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Aorta_Blood_Flow_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Aorta_Blood_Flow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Aorta_Blood_Flow

% Last Modified by GUIDE v2.5 22-May-2018 02:01:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Aorta_Blood_Flow_OpeningFcn, ...
                   'gui_OutputFcn',  @Aorta_Blood_Flow_OutputFcn, ...
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


% --- Executes just before Aorta_Blood_Flow is made visible.
function Aorta_Blood_Flow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Aorta_Blood_Flow (see VARARGIN)

% Choose default command line output for Aorta_Blood_Flow
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
setappdata(0,'handleMainWindow',hObject);


% UIWAIT makes Aorta_Blood_Flow wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Aorta_Blood_Flow_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% Button Exit
function Exit_Callback(hObject, eventdata, handles)
clc;
clear all;




% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrowsePhase1.
% Browse phase images for cross-section 1
function BrowsePhase1_Callback(hObject, eventdata, handles)
% hObject    handle to BrowsePhase1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name  = uigetdir(' ','Select Phase folder 1 (Phase folder)'); %Open folder selection dialog box
files = dir(folder_name); %List folder contents
j=1;
for i = 1:size(files, 1)
    disp(fullfile(folder_name, files(i).name));
    fileName = files(i).name;
    if fileName(1,1) ~= '.'
       handles.imagesPhase1{j} = dicomread(fullfile(folder_name, files(i).name));
       j = j + 1;
    end
end
handles.sliderPhase1.Min = 1;
handles.sliderPhase1.Max = size(handles.imagesPhase1, 2);
handles.sliderPhase1.Value = 1;
handles.sliderPhase1.SliderStep = [1/(size(handles.imagesPhase1, 2)-1) , 10/(size(handles.imagesPhase1, 2)-1)];

axes(handles.axesPhase1)
imshow(handles.imagesPhase1{1}, []);
guidata(hObject, handles);



% --- Executes on button press in BrowseAmplitude1.
% Browse magnitude images for cross-section 1
function BrowseAmplitude1_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------
% This function will open a file dialog in order to display the Magnitude
% images 
%-------------------------------------------------------------
folder_name  = uigetdir(' ','Select magnitude folder 1 (Amplitude folder)'); %Open folder selection dialog box
files = dir(folder_name); %List folder contents
j=1;
for i = 1:size(files, 1)
    disp(fullfile(folder_name, files(i).name));
    fileName = files(i).name;
    if fileName(1,1) ~= '.'
       handles.fileNames{j} = fullfile(folder_name, files(i).name);
       handles.imagesMagnitude1{j} = dicomread(fullfile(folder_name, files(i).name));
       j = j + 1;
    end
end
handles.sliderMagnitude1.Min = 1;
handles.sliderMagnitude1.Max = size(handles.imagesMagnitude1, 2);
handles.sliderMagnitude1.Value = 1;
handles.sliderMagnitude1.SliderStep = [1/(size(handles.imagesMagnitude1, 2)-1) , 10/(size(handles.imagesMagnitude1, 2)-1)];

axes(handles.axesMagnitude1)
imshow(handles.imagesMagnitude1{1}, []);
guidata(hObject, handles);


% --- Executes on button press in selectROI1.
function selectROI1_Callback(hObject, eventdata, handles)
% ----------------------------------------------------------
% This function will select a region of interest, detect a cercle on this
% region 
% Output : center , radius and metric of the cercle with is the ascending
% aorta 

%-----------------------------------------------------------
[xMagnitude1, yMagnitude1] = ginput(1);
if xMagnitude1 > 2 & yMagnitude1 > 1
    handles.pointXMagnitude1 =  xMagnitude1;
    handles.pointYMagnitude1 =  yMagnitude1;
    
    pre_ASC_R=10.25;  % the radius of the first slice ROI ASC circle
    pre_ASC_C=[xMagnitude1 yMagnitude1];  % the center coordinates of the first slice ROI ASC circle
    ASC_R = 10.1;        % the radius of the second slice ROI ASC circle
    ASC_C = [xMagnitude1 yMagnitude1];  % the center coordinates of the second slice ROI ASC circle
    
    ASC_cir = zeros(1,length(handles.imagesMagnitude1));
    ASC_num = zeros(1,length(handles.imagesMagnitude1));
    
    limitASC = 10.7; %distance limit.If the center coordinates distance between the previous slice and the next slice is less than 10,we can accept the error.
    
    for i = 1 : length(handles.imagesMagnitude1)
        [Center, Radius, metric] = imfindcircles(handles.imagesMagnitude1{i},[3 12]);
        [ASC_R,ASC_C,isASC] = segment(Radius,Center,pre_ASC_R,pre_ASC_C,limitASC);
        [ASC_R,ASC_C];
        if isASC==1
            pre_ASC_C(:) = ASC_C(:);
            pre_ASC_R=ASC_R;
        end
        
        figure(1),
        imshow(handles.imagesMagnitude1{i},[]);
        
        if i == get(handles.sliderMagnitude1,'Value')
            cla(handles.axesMagnitude1);
            cla(handles.axesPhase1);
            
            %Show the threshold in Magnitude images
            axes(handles.axesMagnitude1);
            imshow(handles.imagesMagnitude1{i},[])
            viscircles(pre_ASC_C,pre_ASC_R,'EdgeColor','r');
            
            %Show the threshold in Flow images
            axes(handles.axesPhase1);
            imshow(handles.imagesPhase1{i},[])
            viscircles(pre_ASC_C,pre_ASC_R,'EdgeColor','r');
        end
    end
    
    center1 = num2str(ASC_C);
    radius1 = num2str(ASC_R);
set(handles.centerx1, 'String', center1 );
set(handles.radius1, 'String', radius1 );

end

guidata(hObject, handles);





% --- Executes on button press in browsePhase2.
function browsePhase2_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------
% This function will open a file dialog in order to display the Phase
% images 
%-------------------------------------------------------------
folder_name  = uigetdir('','Select phase folder');
files = dir(folder_name);
j = 1;
for i = 1:size(files, 1)
    disp(fullfile(folder_name, files(i).name));
    fileName = files(i).name;
    if fileName(1,1) ~= '.'
        handles.imagesPhase2{j} = dicomread(fullfile(folder_name, files(i).name));
        j = j + 1;
    end
end

handles.sliderPhase2.Min = 1;
handles.sliderPhase2.Max = size(handles.imagesPhase2, 2);
handles.sliderPhase2.Value = 1;
handles.sliderPhase2.SliderStep = [1/(size(handles.imagesPhase2, 2)-1) , 10/(size(handles.imagesPhase2, 2)-1)];

axes(handles.axesPhase2)
imshow(handles.imagesPhase2{1}, [])

guidata(hObject, handles);


% --- Executes on button press in BrowseAmplitude2.
% Browse magnitude images for cross-section 2
function BrowseAmplitude2_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseAmplitude2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name  = uigetdir('','Select magnitude folder');
files = dir(folder_name);
j = 1;
for i = 1:size(files, 1)
    disp(fullfile(folder_name, files(i).name));
    fileName = files(i).name;
    if fileName(1,1) ~= '.'
        handles.imagesMagnitude2{j} = dicomread(fullfile(folder_name, files(i).name));
        j = j + 1;
    end
end

handles.sliderMagnitude2.Min = 1;
handles.sliderMagnitude2.Max = size(handles.imagesMagnitude2, 2);
handles.sliderMagnitude2.Value = 1;
handles.sliderMagnitude2.SliderStep = [1/(size(handles.imagesMagnitude2, 2)-1) , 10/(size(handles.imagesMagnitude2, 2)-1)];

axes(handles.axesMagnitude2)
imshow(handles.imagesMagnitude2{1}, [])

guidata(hObject, handles);

% --- Executes on button press in selectROI2.
function selectROI2_Callback(hObject, eventdata, handles)
% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
[xMagnitude2,yMagnitude2] = ginput(1);
if xMagnitude2 > 1 & yMagnitude2 > 1
    handles.pointXMagnitude2 =  xMagnitude2;
    handles.pointYMagnitude2 =  yMagnitude2;
    
    pre_ASC_R=10.25;  % the radius of the first slice ROI ASC circle
    pre_ASC_C=[xMagnitude2 yMagnitude2];  % the center coordinates of the first slice ROI ASC circle
    ASC_R = 10.1;        % the radius of the second slice ROI ASC circle
    ASC_C = [xMagnitude2 yMagnitude2];  % the center coordinates of the second slice ROI ASC circle
    
    ASC_cir = zeros(1,length(handles.imagesMagnitude2));
    ASC_num = zeros(1,length(handles.imagesMagnitude2));
    
    limitASC = 10.7; %distance limit.If the center coordinates distance between the previous slice and the next slice is less than 10,we can accept the error.
    
    for i = 1 : length(handles.imagesMagnitude2)
        %            img1=imadjust(img0);
        %             figure(1),imshow(handles.imagesMagnitude1{i},[]);
        [Center, Radius, metric] = imfindcircles(handles.imagesMagnitude2{i},[3 12]);
        [ASC_R,ASC_C,isASC] = segment(Radius,Center,pre_ASC_R,pre_ASC_C,limitASC);
        [ASC_R,ASC_C]
        if isASC==1
            pre_ASC_C(:) = ASC_C(:);
            pre_ASC_R=ASC_R;
        end
        
        figure(1),
        imshow(handles.imagesMagnitude2{i},[]);
        
        if i == get(handles.sliderMagnitude2,'Value')
            cla(handles.axesMagnitude2);
            cla(handles.axesPhase2);
            
            %Show the threshold in Magnitude images
            axes(handles.axesMagnitude2);
            imshow(handles.imagesMagnitude2{i},[])
            viscircles(pre_ASC_C,pre_ASC_R,'EdgeColor','r');
            
            %Show the threshold in Flow images
            axes(handles.axesPhase2);
            imshow(handles.imagesPhase2{i},[]);
            viscircles(pre_ASC_C,pre_ASC_R,'EdgeColor','r');
        end
    end
    
end
center2 = num2str(ASC_C);
 radius2 = num2str(ASC_R);
set(handles.center2, 'String', center2 );
set(handles.radius2, 'String', radius2 );
guidata(hObject, handles);


% --- Executes on slider movement.
function sliderPhase1_Callback(hObject, eventdata, handles)
% hObject    handle to sliderPhase1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
axes(handles.axesPhase1);
imshow(handles.imagesPhase1{get(hObject,'Value')}, []);

axes(handles.axesMagnitude1);
magnitudeImage1Handle = imshow(handles.imagesMagnitude1{get(hObject,'Value')}, []);
handles.sliderMagnitude1.Value = get(hObject,'Value');
set(magnitudeImage1Handle,'ButtonDownFcn', @magnitudeImage1ClickCallback);



% --- Executes during object creation, after setting all properties.
function sliderPhase1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderPhase1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderMagnitude1_Callback(hObject, eventdata, handles)
% hObject    handle to sliderMagnitude1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
axes(handles.axesMagnitude1);
imshow(handles.imagesMagnitude1{get(hObject,'Value')}, []);
%    set(magnitudeImage1Handle,'ButtonDownFcn', @magnitudeImage1ClickCallback);

axes(handles.axesPhase1);
imshow(handles.imagesPhase1{get(hObject,'Value')}, []);
handles.sliderPhase1.Value = get(hObject,'Value');


% --- Executes during object creation, after setting all properties.
function sliderMagnitude1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderMagnitude1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderPhase2_Callback(hObject, eventdata, handles)
% hObject    handle to sliderPhase2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
axes(handles.axesPhase2);
imshow(handles.imagesPhase2{get(hObject,'Value')}, []);

axes(handles.axesMagnitude2);
imshow(handles.imagesMagnitude2{get(hObject,'Value')}, []);
handles.sliderMagnitude2.Value = get(hObject,'Value');

% --- Executes during object creation, after setting all properties.
function sliderPhase2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderPhase2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderMagnitude2_Callback(hObject, eventdata, handles)
% hObject    handle to sliderMagnitude2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
axes(handles.axesMagnitude2);
imshow(handles.imagesMagnitude2{get(hObject,'Value')}, []);

axes(handles.axesPhase2);
imshow(handles.imagesPhase2{get(hObject,'Value')}, []);
handles.sliderPhase2.Value = get(hObject,'Value');

% --- Executes during object creation, after setting all properties.
function sliderMagnitude2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderMagnitude2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_Sagital_length.
function pushbutton_Sagital_length_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
% This function estimates the aortic arch length from sagittal plane image
% Output
% length: estimated aortic arch length
%-------------------------------------------------------------------------
% Open a file dialog an choose a file
[filename, pathname] = uigetfile('*.*');
filename = fullfile (pathname, filename);
info = dicominfo(filename);
Y = dicomread(info);
% Display the image
axes(handles.axesSagital);
imshow(Y, 'DisplayRange', []);hold on 

% Initialize the values
coordx = zeros(10,1);
coordy = zeros(10,1);
xs = [];
ys = [];
xold=0;
yold=0;
xc = [];
yc =[];


% Get the measurement given by the user  
MaxMeasures = str2num(get(handles.MeasuresSagital,'String'));

for k = 1:MaxMeasures
    % Draw line between 2 points 
    for h =1:2
        [coordx(h,1), coordy(h,1)] = ginput(1);  
        xs = [xs;coordx(h,1)];
        ys = [ys;coordy(h,1)];
        if xold;
            plot([xold coordx(h,1)],[yold coordy(h,1)],'r.-');
            % Get the center coordinates
            xc = [xc; (xold+coordx(h,1))/2];
            yc = [yc; (yold+coordy(h,1))/2];
        else
            plot(coordx, coordy,'go');
        end 
        xold=coordx(h,1);
        yold=coordy(h,1);
    %coordinates(k,:) = [coordx(k,1), coordy(k,1)];
    %plot(coordinates(:,1), coordinates(:,2), '.-');
    end 
    xold = 0;
end
plot(xc,yc);
ps = info.PixelSpacing(1); % pixel spacing format from DICOM 
length = 0; % initialize the length of the aorta

for k = 1 : MaxMeasures-1
    segment = sqrt(((xc(k,1)-xc(k+1,1))*ps)^2 + ((yc(k,1)-yc(k+1,1))*ps)^2);
    % Show the length of the center points of each lines connected together
    length = length + segment;
    set(handles.SagitalLengthDisplay, 'String', length );
end


% --- Executes on button press in pushbutton_axial_length.
function pushbutton_axial_length_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_axial_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function SagitalLengthDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to SagitalLengthDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SagitalLengthDisplay as text
%        str2double(get(hObject,'String')) returns contents of SagitalLengthDisplay as a double


% --- Executes during object creation, after setting all properties.
function SagitalLengthDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SagitalLengthDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function center1_Callback(hObject, eventdata, handles)
% hObject    handle to center1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of center1 as text
%        str2double(get(hObject,'String')) returns contents of center1 as a double


% --- Executes during object creation, after setting all properties.
function center1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to center1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderSagital_Callback(hObject, eventdata, handles)
% hObject    handle to sliderSagital (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider





% --- Executes during object creation, after setting all properties.
function sliderSagital_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSagital (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function radius1_Callback(hObject, eventdata, handles)
% hObject    handle to radius1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of radius1 as text
%        str2double(get(hObject,'String')) returns contents of radius1 as a double


% --- Executes during object creation, after setting all properties.
function radius1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radius1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MeasuresSagital_Callback(hObject, eventdata, handles)
% hObject    handle to MeasuresSagital (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MeasuresSagital as text
%        str2double(get(hObject,'String')) returns contents of MeasuresSagital as a double


% --- Executes during object creation, after setting all properties.
function MeasuresSagital_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MeasuresSagital (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function center2_Callback(hObject, eventdata, handles)
% hObject    handle to center2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of center2 as text
%        str2double(get(hObject,'String')) returns contents of center2 as a double


% --- Executes during object creation, after setting all properties.
function center2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to center2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function radius2_Callback(hObject, eventdata, handles)
% hObject    handle to radius2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of radius2 as text
%        str2double(get(hObject,'String')) returns contents of radius2 as a double


% --- Executes during object creation, after setting all properties.
function radius2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radius2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function metric2_Callback(hObject, eventdata, handles)
% hObject    handle to metric2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of metric2 as text
%        str2double(get(hObject,'String')) returns contents of metric2 as a double


% --- Executes during object creation, after setting all properties.
function metric2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to metric2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function metric1_Callback(hObject, eventdata, handles)
% hObject    handle to metric1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of metric1 as text
%        str2double(get(hObject,'String')) returns contents of metric1 as a double


% --- Executes during object creation, after setting all properties.
function metric1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to metric1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function centerx1_Callback(hObject, eventdata, handles)
% hObject    handle to centerx1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of centerx1 as text
%        str2double(get(hObject,'String')) returns contents of centerx1 as a double


% --- Executes during object creation, after setting all properties.
function centerx1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to centerx1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function centery_Callback(hObject, eventdata, handles)
% hObject    handle to centery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of centery as text
%        str2double(get(hObject,'String')) returns contents of centery as a double


% --- Executes during object creation, after setting all properties.
function centery_CreateFcn(hObject, eventdata, handles)
% hObject    handle to centery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in segment2.
function segment2_Callback(hObject, eventdata, handles)
% hObject    handle to segment2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pre_ASC_R=10.25;  % the radius of the first slice ROI ASC circle
pre_ASC_C=[handles.pointXMagnitude2 handles.pointYMagnitude2];  % the center coordinates of the first slice ROI ASC circle
ASC_R = 10.1;        % the radius of the second slice ROI ASC circle
ASC_C = [handles.pointXMagnitude2 handles.pointYMagnitude2];  % the center coordinates of the second slice ROI ASC circle

handles.DES_cir = zeros(1,length(handles.imagesMagnitude2));
handles.DES_num = zeros(1,length(handles.imagesMagnitude2));

limitASC = 10.7; %distance limit.If the center coordinates distance between the previous slice and the next slice is less than 10,we can accept the error.

for i = 1 : length(handles.imagesMagnitude1)
    %            img1=imadjust(img0);
    %             figure(1),imshow(handles.imagesMagnitude1{i},[]);
    [Center, Radius, metric] = imfindcircles(handles.imagesMagnitude1{i},[3 12]);
    [ASC_R,ASC_C,isASC] = segment(Radius,Center,pre_ASC_R,pre_ASC_C,limitASC);
    [ASC_R,ASC_C]
    if isASC==1
        pre_ASC_C(:) = ASC_C(:);
        pre_ASC_R=ASC_R;
    end
    
    handles.DesROI{i}=findROI(handles.imagesPhase1{i},pre_ASC_R,pre_ASC_C);
    graylevel=0;
    graylevel=double(graylevel);
    handles.DesROI{i}=double(handles.DesROI{i});
    num=0;
    ROI = handles.DesROI{i};
    for x=1:132
        for y=1:192
            if ROI(x,y)~=0
                graylevel= graylevel+ROI(x,y);
                num= num + 1;
            end
        end
    end
    graylevel=graylevel/num;
    handles.DES_cir(i)= graylevel;
    handles.DES_num(i) = num;
    
    
    
    if i == get(handles.sliderMagnitude2,'Value')
        cla(handles.axesPhase2);
        
        %Show the threshold in Flow images
        axes(handles.axesPhase2);
        imshow(handles.DesROI{i},[])
    end
end

guidata(hObject, handles);


% --- Executes on button press in segment1.
function segment1_Callback(hObject, eventdata, handles)
% hObject    handle to segment1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pre_ASC_R=10.25;  % the radius of the first slice ROI ASC circle
pre_ASC_C=[handles.pointXMagnitude1 handles.pointYMagnitude1];  % the center coordinates of the first slice ROI ASC circle
ASC_R = 10.1;        % the radius of the second slice ROI ASC circle
ASC_C = [handles.pointXMagnitude1 handles.pointYMagnitude1];  % the center coordinates of the second slice ROI ASC circle

handles.ASC_cir = zeros(1,length(handles.imagesMagnitude1));
handles.ASC_num = zeros(1,length(handles.imagesMagnitude1));

limitASC = 10.7; %distance limit.If the center coordinates distance between the previous slice and the next slice is less than 10,we can accept the error.

for i = 1 : length(handles.imagesMagnitude1)
    [Center, Radius, metric] = imfindcircles(handles.imagesMagnitude1{i},[3 12]);
    [ASC_R,ASC_C,isASC] = segment(Radius,Center,pre_ASC_R,pre_ASC_C,limitASC);
    [ASC_R,ASC_C]
    if isASC==1
        pre_ASC_C(:) = ASC_C(:);
        pre_ASC_R=ASC_R;
    end
    
    handles.AscROI{i}=findROI(handles.imagesPhase1{i},pre_ASC_R,pre_ASC_C);
    graylevel=0;
    graylevel=double(graylevel);
    handles.AscROI{i}=double(handles.AscROI{i});
    num=0;
    ROI = handles.AscROI{i};
    for x=1:132
        for y=1:192
            if ROI(x,y)~=0
                graylevel= graylevel+ROI(x,y);
                num= num + 1;
            end
        end
    end
    graylevel=graylevel/num;
    handles.ASC_cir(i)= graylevel;
    handles.ASC_num(i) = num;
    
    
    
    if i == get(handles.sliderMagnitude1,'Value')
        cla(handles.axesPhase1);
        
        %Show the threshold in Flow images
        axes(handles.axesPhase1);
        imshow(handles.AscROI{i},[])
    end
end

guidata(hObject, handles);


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Flow.
function Flow_Callback(hObject, eventdata, handles)
% hObject    handle to Flow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Flow_1 = zeros(1,length(handles.imagesPhase1));
handles.Flow_2 = zeros(1,length(handles.imagesPhase2));

for i=1:length(handles.imagesPhase1)
    handles.Flow_1(i) = handles.Speed_1(i)/100 * 1.67^2 * handles.ASC_num(i); %Flow=speed * area
    handles.Flow_2(i) = handles.Speed_2(i)/100 * 1.67^2 * handles.DES_num(i);
end

time = zeros(1, length(handles.imagesPhase1));

for i=1:length(handles.fileNames)
    info{i}=dicominfo(handles.fileNames{i});
    time(i) = info{i}.TriggerTime; 
end

assignin('base', 'time', time);
axes(handles.axesFlow);
plot(time,handles.Flow_1,'r'),hold on;
plot(time,handles.Flow_2,'b');
title('Blood Flow Volume(ml/sec)');
legend('Ascending aorta','Descending aorta','Best')
xlabel('Trigger delay (ms)')
hold off;
guidata(hObject, handles);


% --- Executes on button press in Velocity.
function Velocity_Callback(hObject, eventdata, handles)
% ----------------------------------------------------------
% This function will return the velocity(speed)
%
handles.Speed_1=zeros(1,length(handles.imagesPhase1));
handles.Speed_2=zeros(1,length(handles.imagesPhase2));

NV = 2048 ;
MV = 150 ;

for i=1:length(handles.imagesPhase1)
    SI_asc = handles.ASC_cir(i);
    SI_dsc = handles.DES_cir(i);
    % Equation [2] for the descending aorta.
    handles.Speed_1(i) = ((SI_asc - NV)*MV) / NV; %Calculate the pixels and get the speed
    % Equation[1] for the ascending aorta
    handles.Speed_2(i) = ((NV - SI_dsc ) * MV) / NV;
end

time = zeros(1, length(handles.imagesPhase1));

for i=1:length(handles.fileNames)
    info{i}=dicominfo(handles.fileNames{i});
    time(i) = info{i}.TriggerTime; 
end

assignin('base', 'time', time);

axes(handles.axesVelocity);
plot(time,handles.Speed_1,'r'),hold on;
plot(time,handles.Speed_2,'b');
title('Velocity(cm/sec)');
legend('Ascending aorta','Descending aorta','Best')
xlabel('Trigger delay (ms)')
hold off;
guidata(hObject, handles);
