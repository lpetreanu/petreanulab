function varargout = Leo_Ca_Analyzer_SI2019a_EP(varargin)
% LEO_CA_ANALYZER M-file for Leo_Ca_Analyzer_SI2019a_EP.fig
%      Leo_Ca_Analyzer_SI2019a_EP, by itself, creates a new LEO_CA_ANALYZER or raises the existing
%      singleton*.
%
%      H = Leo_Ca_Analyzer_SI2019a_EP returns the handle to a new LEO_CA_ANALYZER or the handle to
%      the existing singleton*.
%
%      Leo_Ca_Analyzer_SI2019a_EP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Leo_Ca_Analyzer_SI2019a_EP.M with the given input arguments.
%
%      Leo_Ca_Analyzer_SI2019a_EP('Property','Value',...) creates a new Leo_Ca_Analyzer_SI2019a_EP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Leo_Ca_Analyzer_SI2019a_EP_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Leo_Ca_Analyzer_SI2019a_EP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% *************************************************************************
%
%Bug Fix workaround for 2007b
%
% This bug was fixed as of R2008a.
%
% If you have a current subscription to MathWorks Software Maintenance Service (SMS), you can download product updates. If not, learn more about MathWorks SMS.
% Workaround
% If you are using a previous version, please read the following:
%
% To work around this issue, replace the following files in your existing MATLAB installation with the updated versions attached.
%
%    1. Determine the MATLAB root directory on your machine by typing
%
%       matlabroot
%
%       at the MATLAB prompt. This will be referred to as $MATLAB below.
%    2. Quit MATLAB.
%    3. Move to the following directory:
%
%       cd $MATLAB/toolbox/images/imuitools/private
%    4. Make a back up copy of the following files.
%
%       mv polygonSymbol.m polygonSymbol.m.old
%       mv manageInteractivePlacement.m manageInteractivePlacement.m.old
%
%    5. Download the attached M-files and store them in the current directory.
%    6. Restart MATLAB
%**************************************************************************
%
% Edit the above text to modify the response to help Leo_Ca_Analyzer_SI2019a_EP

% Last Modified by GUIDE v2.5 13-Aug-2020 13:41:46

%% EDITS BY TIAGO
%   1. Removed auto-advance ROI number from the setDonutROI function


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Leo_Ca_Analyzer_SI2019a_EP_OpeningFcn, ...
    'gui_OutputFcn',  @Leo_Ca_Analyzer_SI2019a_EP_OutputFcn, ...
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


% --- Executes just before Leo_Ca_Analyzer_SI2019a_EP is made visible.

function Leo_Ca_Analyzer_SI2019a_EP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Leo_Ca_Analyzer (see VARARGIN)

% Choose default command line output for Leo_Ca_Analyzer
handles.output = hObject;
%handles.timePerLine = 0.002; % in s
handles.roiNumber=1;
colist=rand(20,3)*0.9+0.1;
handles.colorlist=repmat(colist,[50 1]);
handles.displayChannel='G';
handles.imageHandle=[];
handles.currentColorMapMaxValue=500;
handles.currentColorMapMinValue=0;
handles.analysisMode='F';
handles.planeOfInterest=1;
handles.bidi_correct_offset_flag = 0;
handles.bidiOffset = 0;
%set(handles.axes1,'ButtonDownFcn',{@ImClickCallback, handles});
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Leo_Ca_Analyzer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Leo_Ca_Analyzer_SI2019a_EP_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
try
    delete(handles.imageHandle)
catch
end
availablePlanes=get(handles.PlaneNumber,'String');
if handles.state.fastZEnable
    handles.planeOfInterest=str2double(availablePlanes(get(handles.PlaneNumber,'Value')));
else
    handles.planeOfInterest=1;
end


switch handles.displayChannel
    case 'G'
        if handles.numberOfPlanes>1
            activeChannel=handles.channels.G{handles.planeOfInterest};
        else
            activeChannel=handles.channels.G;
        end
    case 'R'
        if handles.numberOfPlanes>1
            activeChannel=handles.channels.R{handles.planeOfInterest};
        else
            activeChannel=handles.channels.R;
        end
    case 'ProjectAcrossTrials'
        if ~isfield(handles,'trialsProjection')
            disp('Calculate Projection across trials first. Switching to Green Single Frame')
            if handles.numberOfPlanes>1
                activeChannel=handles.channels.G{handles.planeOfInterest};
            else
                activeChannel=handles.channels.G;
            end
        else
            activeChannel=handles.trialsProjection;
        end
end

frame=round(get(hObject,'Value'));
im = showImage(activeChannel(:,:,frame),handles);

handles.imageHandle=im;
colormap(gray);
%set(handles.imageHandle,'ButtonDownFcn',{@ImClickCallback, handles});
uistack(handles.imageHandle,'bottom')
frameTime=handles.numberOfRows * handles.timePerLine*frame*handles.numberOfPlanes;
set(handles.FrameCounter,'String',[num2str(frame) '/' num2str(handles.numberOfFrames) ' | '  sprintf('%0.3f',frameTime) 's']);
guidata(handles.figure1,handles) ;
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in loadButton.
function loadButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'imageHandle')
    delete(handles.imageHandle)
    
end
if isfield(handles,'channels')
    handles.channels=[];
end
if isfield(handles,'ccimage');
    handles=rmfield(handles, 'ccimage');
end
[f,p]=uigetfile('*.tif','Choose file to load');
cd(p)
handles.currentFile=f;
handles.currentDirectory=p;



try
    [H,Data] = load_SI2019a_raw(p,f);
    handles.state=H;
    handles.timePerLine=H.scanLinePeriod;
    numPlanes=H.stackNumSlices;
    numChannels=size(H.channelsActive,1);
    Data(Data==0)=NaN;
    handles.pmt_offset=H.SI.hChannels.channelOffset;
    
    if numPlanes>1
        if numChannels>1
            for i=1:numPlanes;
%                 handles.channels.G{i}=squeeze(Data(:,:,1,:,i,i:numPlanes:end));
%                 handles.channels.R{i}=squeeze(Data(:,:,2,:,i,i:numPlanes:end));
                handles.channels.G{i}=squeeze(Data(:,:,1,:,i,:));
                handles.channels.R{i}=squeeze(Data(:,:,2,:,i,:));
            end
        elseif numChannels<2
            for i=1:numPlanes;
%                 handles.channels.G{i}=squeeze(Data(:,:,1,:,i,i:numPlanes:end));
                handles.channels.G{i}=squeeze(Data(:,:,1,:,i,:));
            end
        end
    else
        handles.channels.G=squeeze(Data(:,:,1,:));
    end
    
catch % in case of opening files saved with Efthychios' registration (only works for one plane, one channel)
   [H,Data] = load_SI2019a_EP(p,f);

    handles.timePerLine=H.scanLinePeriod;

    numPlanes=1;
    Data(Data==0)=NaN;
    handles.channels.G=Data;
    handles.state=H;

end

handles.numberOfPlanes = numPlanes;

if get(handles.loadAvgFilterCheck,'Value')==1
    filter=fspecial('average',[2 2]);
    disp('filtering....')
    fImage=imfilter(handles.channels.G, filter, 'conv');
    handles.channels.G=fImage;
end
% s=size(Data);
% if s(3)==2
%     handles.channels.R=squeeze(Data(:,:,2,:));
%     if get(handles.loadAvgFilterCheck,'Value')==1
%         filter=fspecial('average',[2 2]);
%         disp('filtering....')
%         fImage=imfilter(handles.channels.R, filter, 'conv');
%         handles.channels.R=fImage;
%     end
% end

if iscell(handles.channels.G)
    [r, c z]=size(handles.channels.G{i});
else
    [r, c z]=size(handles.channels.G);
end
handles.numberOfFrames=z;
handles.numberOfRows=r;
handles.numberOfColumns=c;
set(handles.figure1,'Name',handles.currentFile);
set(handles.slider1,'Min',1,'Max',z,'Value',1,'SliderStep',[ 1/(z-1) 10/(z-1)]);
set(handles.FrameCounter,'String',[num2str(1) '/' num2str(z)]);
set(handles.zoomFactorBox,'String', ['Zoom: ' num2str(handles.state.scanZoomFactor)]);
if ~isnan(handles.state.date)
set(handles.DateBox,'String', ['Date: '  handles.state.date]);
end
set(handles.PowerBox,'String', ['Power: ' num2str(handles.state.beamPowers) '%']);
set(handles.FrequencyBox,'String', ['Acq Freq: ' num2str(handles.state.scanFrameRate) ' Hz']);
if handles.state.fastZEnable
    set(handles.IsVolumetricImagingBox,'String',['Volumetric Imaging ON']);
    set(handles.DistanceAcrossPlanesBox,'String',['InterPlane Dist: ', num2str(handles.state.stackZStepSize) ' um']);
    
else
    set(handles.IsVolumetricImagingBox,'String',['Volumetric Imaging OFF']);
    set(handles.DistanceAcrossPlanesBox,'String',['InterPlane Dist: ---']);
    handles.planeOfInterest=1;
end
set(handles.NumberOfImagedPlanesBox,'String',['Number of Planes: ' num2str(handles.state.stackNumSlices) ]);
set(handles.PlaneNumber,'String',[1:handles.state.stackNumSlices],'Value',handles.planeOfInterest);
% handles.timePerLine=handles.state.scanLinePeriod;
guidata(handles.figure1,handles) ;

% --- Executes on button press in BGButton.
function BGButton_Callback(hObject, eventdata, handles)
% hObject    handle to BGButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%[bw, xbw, ybw]=roipolyold;
zoom off
[bw, xbw, ybw]=roipoly;
try
    delete(handles.backgroundHandle)
catch
end

hold on
bghandle=plot(xbw, ybw,'w-');
set(bghandle,'Tag','BG');

handles.backgroundpolyCoord=[ xbw ybw];
handles.backgroundHandle=bghandle;
guidata(handles.figure1,handles) ;

% --- Executes on button press in ROIButton.
function ROIButton_Callback(hObject, eventdata, handles)
% hObject    handle to ROIButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.uitoggletool1,'State','off');
set(handles.uitoggletool3,'State','off');
zoom off
% zoom reset

roiActive=get(handles.listboxROIs,'Value');
colorList=get(handles.axes1,'ColorOrder');
[roi roixActive roiyActive]=roipoly;

try
    delete(handles.roiPlotHandle{roiActive},handles.roiTextHandle(roiActive),handles.ROIPixelList{roiActive}.nuclear )
catch
end

%roiActive=repmat(roiActive, [ 1 1 handles.numberOfFrames]);
%roiActive=uint8(roi);
hold on

handles.roiPlotHandle{roiActive}=plot(roixActive, roiyActive,'Color', handles.colorlist(roiActive,:), 'LineWidth',2);
%set(handles.roiPlotHandle(roiActive),'Tag','ROI','ButtonDownFcn','Leo_moverotateobj(''BtnDown'')');
htext=text(mean(roixActive), mean(roiyActive), num2str(roiActive),'Color',handles.colorlist(roiActive,:),'FontWeight','bold','FontSize',6);
handles.roiTextHandle(roiActive)=htext;
handles.roipolyCoord{roiActive}.outer=[ roixActive roiyActive];
handles.roipolyCoord{roiActive}.inner=[]; % not implemented yet
handles.ROIPixelList{roiActive}.outer=find(roi);
handles.ROIPixelList{roiActive}.inner=[];

guidata(handles.figure1,handles) ;

% --- Executes on button press in AnalyzeButton.
function AnalyzeButton_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%roiActive=get(handles.listboxROIs,'Value')
tic
switch handles.displayChannel
    case 'G'
        if handles.numberOfPlanes>1
            activeChannel=handles.channels.G{handles.planeOfInterest};
        else
            activeChannel=handles.channels.G;
        end
    case 'R'
        if handles.numberOfPlanes>1
            activeChannel=handles.channels.R{handles.planeOfInterest};
        else
            activeChannel=handles.channels.R;
        end
    case 'ProjectAcrossTrials'
        fprintf('ATTN! Please select Green or Red as the display channels before analyzing\n\n\n')
end

% if bidirection correction is on, analyze the corrected frames
if handles.bidi_correct_offset_flag && handles.bidiOffset ~= 0
    activeChannel=biDiOffsetCorrection(activeChannel,handles.bidiOffset);
end

bw=poly2mask(handles.backgroundpolyCoord(:,1),handles.backgroundpolyCoord(:,2),handles.numberOfRows,handles.numberOfColumns);
pixNumRoiBG=sum(bw(:));
bw=repmat(bw, [ 1 1 handles.numberOfFrames]);
if strcmp(class(activeChannel),'uint8')
    bw=uint8(bw);
elseif strcmp(class(activeChannel),'int16')
    bw=int16(bw);
else
    bw=uint16(bw);
end
handles.backgroundpoly=bw;


background=handles.backgroundpoly.*activeChannel;
%background=double(background);
%background(background==0)=NaN;
areaOfInterestBG=reshape(activeChannel(logical(bw)),pixNumRoiBG, handles.numberOfFrames);
%areaOfInterestBG=double(areaOfInterestBG);
backgroundTrace=mean(areaOfInterestBG,1);
B=nanmean(background(logical(bw)));

handles.background=B;
handles.lastBackgroundTrace=backgroundTrace;

contourWidth=str2double(get(handles.donutCellRadius,'String'));
se90 = strel('line', contourWidth/2, 90); %extended neuron border
se0 = strel('line', contourWidth/2, 0);

se90N = strel('line', contourWidth*30, 90); % Neuropil radios
se0N = strel('line', contourWidth*30, 0);
%% Project all ROI to calculate Neuropil later
if ~isfield(handles,'currentTrialNumber');
    handles.startTrialNum=NaN;
    handles.currentTrialNumber=NaN;
end

if handles.currentTrialNumber==handles.startTrialNum||hObject==handles.AnalyzeButton;
    if get(handles.CalculateNeuropil,'Value')
        %FFContour=zeros(handles.roiNumber,handles.numberOfFrames/handles.numberOfPlanes);
        ROIContourMaskAllROIs=false(handles.numberOfRows,handles.numberOfColumns,handles.roiNumber);
        
        
        AllRoisMask= false(handles.numberOfRows,handles.numberOfColumns);
        AllGroiMask= false(handles.numberOfRows,handles.numberOfColumns,handles.roiNumber);
        
        AllNuclearRoiMask= false(handles.numberOfRows,handles.numberOfColumns,handles.roiNumber);
        for i=1:handles.roiNumber; % Calculate All ROIs Surfaces for Neuropil
            
            GroiMask=false(handles.numberOfRows,handles.numberOfColumns);
            NuclearRoiMask=false(handles.numberOfRows,handles.numberOfColumns);
            GroiMask(handles.ROIPixelList{i}.outer)=1;
            if isfield(handles.ROIPixelList{i},'nuclear')
                NuclearRoiMask(handles.ROIPixelList{i}.nuclear)=1;
            end
            AllGroiMask(:,:,i)=GroiMask;
            AllNuclearRoiMask(:,:,i)=NuclearRoiMask;
            AllRoisMask=AllRoisMask|GroiMask+NuclearRoiMask;
            
        end
        AllRoisMaskDilated=imdilate(AllRoisMask,[ se90 se0]);
        handles.AllRoisMaskDilated=AllRoisMaskDilated;
    end
end
%%

if get(handles.UseGPU,'Value')
    
    
    % AllNuclearRoiMaskCPU=gather(AllRoisMask);
    %activeChannelGPU=gpuArray(activeChannel);
    
    
    FF=gpuArray.zeros(handles.roiNumber,handles.numberOfFrames/handles.numberOfPlanes);
    for i=1:handles.roiNumber;
        
        %GroiMaskGPU=gpuArray.zeros(handles.numberOfRows,handles.numberOfColumns);
        if isempty(handles.ROIPixelList{i}.outer)  %patch for bug where saved ROIs have missing pixel list
            handles.ROIPixelList{i}.outer=poly2mask(handles.roipolyCoord{i}.outer(:,1),handles.roipolyCoord{i}.outer(:,2),handles.numberOfRows,handles.numberOfColumns);
        end
        
        % GPU
        
        %         roiMask=gpuArray((AllGroiMask(:,:,i))); %%Calculate F
        %         pixNumRoi=sum(roiMask(:));
        %        %GroiMaskGPU=logical(repmat(roiMask, [ 1 1 handles.numberOfFrames/handles.numberOfPlanes]));
        %         GroiMaskGPU=repmat(roiMask, [ 1 1 handles.numberOfFrames/handles.numberOfPlanes]);
        %        areaOfInterest=reshape(activeChannel(GroiMaskGPU),pixNumRoi, handles.numberOfFrames/handles.numberOfPlanes);
        %         FF(i,:)=mean(areaOfInterest,1);
        %
        %         if get(handles.CalculateNeuropil,'Value')  % Calculate neuropil
        %         roiMaskcpu=gather(roiMask);
        %         ROIdil=imdilate(roiMaskcpu,[ se90 se0] );
        %         ROIContourMask=(logical(ROIdil-roiMaskcpu))&~(AllRoisMask);
        %         ROIContourMaskAllROIs(:,:,i)=ROIContourMask;
        %         ROIContourMaskGPU=gpuArray(ROIContourMask);
        %         pixNumRoiCont=sum(ROIContourMaskGPU(:));
        %         ContourRoiMask=logical(repmat(ROIContourMaskGPU, [ 1 1 handles.numberOfFrames/handles.numberOfPlanes]));
        %         areaOfInterest=reshape(activeChannel(ContourRoiMask),pixNumRoiCont, handles.numberOfFrames/handles.numberOfPlanes);
        %         FFContour(i,:)=mean(areaOfInterest,1);
        %         FNeuropil=gather(FFContour);
        %         end
        
        
        
        %% newMethod Without Big Repmat
        
        roiMaskOneFrame=gpuArray.zeros(handles.numberOfRows,handles.numberOfColumns);
        
        roiMaskOneFrame(handles.ROIPixelList{i}.outer)=1;
        pixNumRoi=sum(roiMaskOneFrame(:));
        pixInd=find(roiMaskOneFrame);
        totalPixPerFrame=numel(roiMaskOneFrame);
        %Pix Ind of al lframes based on the size of each Frame
        pixIndRep=repmat([pixInd],[1 handles.numberOfFrames/handles.numberOfPlanes]);
        indAllFrames=repmat(0:totalPixPerFrame:totalPixPerFrame*(handles.numberOfFrames/handles.numberOfPlanes-1),[pixNumRoi 1]);
        roiMaskAllFramesInd=pixIndRep+indAllFrames;
        
        areaOfInterest=reshape(activeChannel(roiMaskAllFramesInd(:)),pixNumRoi, handles.numberOfFrames/handles.numberOfPlanes);
        FF(i,:)=mean(areaOfInterest,1);
        
        
        %% Neuropil
        
        if get(handles.CalculateNeuropil,'Value')  % Calculate neuropil
            roiMaskOneFrame=AllGroiMask(:,:,i);
            ROIdil=imdilate(roiMaskOneFrame,[ se90 se0] );
            ROIContourMask=(logical(ROIdil-roiMaskOneFrame))&~(AllRoisMask);
            ROIContourMaskAllROIs(:,:,i)=ROIContourMask;
            
            
            pixNumRoiCont=sum(ROIContourMask(:));
            roiContInd=find(ROIContourMask);
            roiContPixIndRep=repmat([roiContInd],[1 handles.numberOfFrames/handles.numberOfPlanes]);
            roiContIndAllFrames=repmat(0:totalPixPerFrame:totalPixPerFrame*(handles.numberOfFrames/handles.numberOfPlanes-1),[pixNumRoiCont 1]);
            ContourRoiMaskInd=roiContPixIndRep+roiContIndAllFrames;
            
            
            %ContourRoiMask=logical(repmat(ROIContourMask, [ 1 1 handles.numberOfFrames/handles.numberOfPlanes]));
            areaOfInterestN=reshape(activeChannel(ContourRoiMaskInd(:)),pixNumRoiCont, handles.numberOfFrames/handles.numberOfPlanes);
            FNeuropil(i,:)=mean(areaOfInterestN,1);
        end
        
    end
    F=gather(FF);
    
    
else  %CPU
    F=zeros(handles.roiNumber,handles.numberOfFrames/handles.numberOfPlanes);
    for i=1:handles.roiNumber;
        
        %         %     if isempty(handles.ROIPixelList{i}.outer)  %patch for bug where saved ROIs have missing pixel list
        %         %         handles.ROIPixelList{i}.outer=poly2mask(handles.roipolyCoord{i}.outer(:,1),handles.roipolyCoord{i}.outer(:,2),handles.numberOfRows,handles.numberOfColumns);
        %         %     end
        %         roiMaskOneFrame=zeros(handles.numberOfRows,handles.numberOfColumns);
        %         roiMaskOneFrame(handles.ROIPixelList{i}.outer)=1;
        %         %roiMask=poly2mask(handles.roipolyCoord{i}(:,1),handles.roipolyCoord{i}(:,2),handles.numberOfRows,handles.numberOfColumns);
        %         pixNumRoi=sum(roiMaskOneFrame(:));
        %         roiMaskAllFrames=logical(repmat(roiMaskOneFrame, [ 1 1 handles.numberOfFrames/handles.numberOfPlanes]));
        %         areaOfInterest=reshape(activeChannel(roiMaskAllFrames),pixNumRoi, handles.numberOfFrames/handles.numberOfPlanes);
        %         F(i,:)=mean(areaOfInterest,1);
        
        %% newMethod Without Big Repmat
        
        roiMaskOneFrame=zeros(handles.numberOfRows,handles.numberOfColumns);
        roiMaskOneFrame(handles.ROIPixelList{i}.outer)=1;
        pixNumRoi=sum(roiMaskOneFrame(:));
        pixInd=find(roiMaskOneFrame);
        totalPixPerFrame=numel(roiMaskOneFrame);
        %Pix Ind of al lframes based on the size of each Frame
        pixIndRep=repmat([pixInd],[1 handles.numberOfFrames/handles.numberOfPlanes]);
        indAllFrames=repmat(0:totalPixPerFrame:totalPixPerFrame*(handles.numberOfFrames/handles.numberOfPlanes-1),[pixNumRoi 1]);
        roiMaskAllFramesInd=pixIndRep+indAllFrames;
        
        areaOfInterest=reshape(activeChannel(roiMaskAllFramesInd(:)),pixNumRoi, handles.numberOfFrames/handles.numberOfPlanes);
        F(i,:)=mean(areaOfInterest,1);
        
        
        %% Neuropil
        
        if get(handles.CalculateNeuropil,'Value')  % Calculate neuropil
            
            % ROIdil=imdilate(roiMaskOneFrame,[ se90 se0] );
            ROIdilN=imdilate(roiMaskOneFrame,[ se90N se0N] );
            ROIContourMask=(logical(ROIdilN&~handles.AllRoisMaskDilated));
            ROIContourMaskAllROIs(:,:,i)=ROIContourMask;
            
            
            pixNumRoiCont=sum(ROIContourMask(:));
            roiContInd=find(ROIContourMask);
            roiContPixIndRep=repmat([roiContInd],[1 handles.numberOfFrames/handles.numberOfPlanes]);
            roiContIndAllFrames=repmat(0:totalPixPerFrame:totalPixPerFrame*(handles.numberOfFrames/handles.numberOfPlanes-1),[pixNumRoiCont 1]);
            ContourRoiMaskInd=roiContPixIndRep+roiContIndAllFrames;
            
            
            %ContourRoiMask=logical(repmat(ROIContourMask, [ 1 1 handles.numberOfFrames/handles.numberOfPlanes]));
            areaOfInterestN=reshape(activeChannel(ContourRoiMaskInd(:)),pixNumRoiCont, handles.numberOfFrames/handles.numberOfPlanes);
            FNeuropil(i,:)=mean(areaOfInterestN,1);
            
        end
        
    end
end


timePerLine = handles.timePerLine; % in s
deltaT = handles.numberOfRows * timePerLine*handles.numberOfPlanes; % in s, time per frame

switch handles.analysisMode
    case 'F'
        %F=F;%;-B;
        yFactor=50;
        offsetVector=yFactor*(0:handles.roiNumber-1)';
        offsetArray=repmat(offsetVector, [1 handles.numberOfFrames/handles.numberOfPlanes ]);
        Fplot=F+offsetArray;
        
        YLabel='F';
        
        
    case 'deltaF/F'
        %B=0;
        Fo=(prctile((F-B)',30))'; % using lower 30 percentile as Fo!
        
        %Fo=median(F-B,2);
        FoArray=repmat(Fo,[1 handles.numberOfFrames/handles.numberOfPlanes]);
        F=(F-B-FoArray)./FoArray;
        
        yFactor=1;    %<-----Offset between traces for the plot
        offsetVector=yFactor*(0:handles.roiNumber-1)';
        offsetArray=repmat(offsetVector, [1 handles.numberOfFrames/handles.numberOfPlanes ]);
        Fplot=F+offsetArray;
        YLabel='deltaF/F';
    case 'deltaF/F/R'
        
        
        
        backgroundG=handles.backgroundpoly.*handles.channels.G;
        backgroundR=handles.backgroundpoly.*handles.channels.R;
        backgroundG=double(backgroundG);
        backgroundR=double(backgroundR);
        %background(background==0)=NaN;
        
        
        BG=nanmean(backgroundG(logical(bw)));
        BR=nanmean(backgroundR(logical(bw)));
        
        F=zeros(handles.roiNumber,handles.numberOfFrames);
        FR=zeros(handles.roiNumber,handles.numberOfFrames);
        
        for i=1:handles.roiNumber;
            
            
            roiMask=poly2mask(handles.roipolyCoord{i}(:,1),handles.roipolyCoord{i}(:,2),handles.numberOfRows,handles.numberOfColumns);
            pixNumRoi=sum(roiMask(:));
            roiMask=repmat(roiMask, [ 1 1 handles.numberOfFrames]);
            
            areaOfInterestG=reshape(handles.channels.G(roiMask),pixNumRoi, handles.numberOfFrames);
            areaOfInterestR=reshape(handles.channels.R(roiMask),pixNumRoi, handles.numberOfFrames);
            areaOfInterestG=double(areaOfInterestG);
            areaOfInterestR=double(areaOfInterestR);
            
            F(i,:)=mean(areaOfInterestG,1);
            FR(i,:)=mean(areaOfInterestR,1);
            
            
            
        end
        
        
        Red=(FR-BR);
        
        Fo=(prctile(((F-BG)./Red)',15))'; % using lower 35 percentile as Fo!
        FoArray=repmat(Fo,[1 handles.numberOfFrames]);
        F=((F-BG-FoArray)./FoArray)./Red;
        
        yFactor=1;    %<-----Offset between traces for the plot
        offsetVector=yFactor*(0:handles.roiNumber-1)';
        offsetArray=repmat(offsetVector, [1 handles.numberOfFrames ]);
        Fplot=F+offsetArray;
        YLabel='deltaF/F/R';
end

if hObject==handles.AnalyzeButton
    h=figure('Position',[260         0         500        1000]);
    set(h,'DefaultAxesColorOrder',handles.colorlist);
    hax=axes;
    
    set(gcf,'Name', handles.currentFile)
    set(h,'Name', handles.currentFile);
    
    plot(deltaT:deltaT:(handles.numberOfFrames/handles.numberOfPlanes)*deltaT, Fplot);
    
    set(get(hax,'XLabel'),'String', 'Time (s)');
    set(get(hax,'YLabel'),'String', YLabel);
    set(hax,'XLim',[0 handles.numberOfFrames*deltaT])%, 'YLim', [-.1 handles.roiNumber+2]);
    assignin('base', ['F_' handles.currentFile(1:end-4)], F);
    if get(handles.CalculateNeuropil,'Value')
        assignin('base', ['F_Neuropil_' handles.currentFile(1:end-4)], FNeuropil);
        handles.ROIContourMaskAllROIs=ROIContourMaskAllROIs;
        guidata(handles.figure1,handles) ;
    end
else
    handles.lastF=zeros(handles.roiNumber,handles.numberOfFrames);
    handles.lastF=F;
    if get(handles.CalculateNeuropil,'Value')
        handles.lastFNeuropil=FNeuropil;
        handles.ROIContourMaskAllROIs=ROIContourMaskAllROIs;
    else
        handles.lastFNeuropil=NaN;
    end
    %handles.lastF(1,1)
    guidata(handles.figure1,handles) ;
    
end
toc

%

% --- Executes on slider movement.
function ColormapMax_Callback(hObject, eventdata, handles)
% hObject    handle to ColormapMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

maxImageValue=round(get(hObject,'Value'));
set(gca,'CLim', [handles.currentColorMapMinValue maxImageValue ])
set(handles.currentColorMapMax,'String',num2str(maxImageValue));
handles.currentColorMapMaxValue=maxImageValue;
guidata(handles.figure1,handles) ;

% --- Executes during object creation, after setting all properties.
function ColormapMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ColormapMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function MinValueText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinValueText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on slider movement.
function ColormapMin_Callback(hObject, eventdata, handles)
% hObject    handle to ColormapMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
minImageValue=round(get(hObject,'Value'));
set(gca,'CLim', [minImageValue handles.currentColorMapMaxValue]);
set(handles.currentColorMapMin,'String', num2str(minImageValue));
handles.currentColorMapMinValue= minImageValue;
guidata(handles.figure1,handles) ;

% --- Executes during object creation, after setting all properties.
function ColormapMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ColormapMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on button press in MaxProjection.
function MaxProjection_Callback(hObject, eventdata, handles)


% hObject    handle to MaxProjection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.figure1,'Name',handles.currentFile);
% set(handles.ColormapMax,'Min',0, 'Max',max(max(max(handles.channels.G,[],3))), 'Value', 200,...
%     'SliderStep', [ 1/double(max(handles.channels.G(:))) 10/double(max(handles.channels.G(:)))])%
% set(handles.ColormapMin,'Min',0, 'Max',max(max(max(handles.channels.G,[],3))), 'Value',0,...
%     'SliderStep', [ 1/double(max(handles.channels.G(:))) 10/double(max(handles.channels.G(:)))])%
% set(handles.MaxValueText, 'String', [ '/' num2str(max(max(max(handles.channels.G,[],3))))]);
% set(handles.MinValueText, 'String', ['/' num2str(min(min(min(handles.channels.G,[],3))))]);
try
    delete(handles.imageHandle)
    
catch
end


switch handles.displayChannel
    case 'G'
        if handles.numberOfPlanes>1
            activeChannel=handles.channels.G{handles.planeOfInterest};
        else
            activeChannel=handles.channels.G;
        end
    case 'R'
        if handles.numberOfPlanes>1
            activeChannel=handles.channels.R{handles.planeOfInterest};
        else
            activeChannel=handles.channels.R;
        end
    case 'ProjectAcrossTrials'
        if ~isfield(handles,'trialsProjection')
            disp('Calculate Projection across trials first. Switching to Green Single Frame')
            if handles.numberOfPlanes>1
                activeChannel=handles.channels.G{handles.planeOfInterest};
            else
                activeChannel=handles.channels.G;
            end
        else
            activeChannel=handles.trialsProjection;
        end
end



movAvgim=im_mov_avg(activeChannel,5);%
maxZ=max(movAvgim,[],3);
maxValue=double(max(maxZ(:)));
(handles.axes1);


im = showImage(maxZ,handles);
colormap(gray);
handles.imageHandle=im;
%set(handles.imageHandle,'ButtonDownFcn',{@ImClickCallback, handles});
set(handles.ColormapMax,'Min',0, 'Max',maxValue, 'Value', maxValue/2,'SliderStep',[ 1/maxValue 10/maxValue]);%
set(handles.ColormapMin,'Min',min(min(maxZ(:))), 'Max',maxValue, 'Value',min(min(maxZ(:))), 'SliderStep', [ 1/maxValue 10/maxValue]);
set(handles.MaxValueText, 'String', [ '/' num2str(maxValue)]);
set(handles.MinValueText, 'String', ['/' num2str(min(maxZ(:)))]);
set(handles.currentColorMapMax,'String',num2str(maxValue/2));


uistack(handles.imageHandle,'bottom')
guidata(handles.figure1,handles) ;



% --- Executes on button press in CopyImageButton.
function CopyImageButton_Callback(hObject, eventdata, handles)
% hObject    handle to CopyImageButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
axes
set(gca,'YDir','reverse')
try
    copyHandle=copyobj([handles.imageHandle handles.roiPlotHandle handles.roiTextHandle],gca);
catch
    copyHandle=copyobj([handles.imageHandle],gca);
end
set(gca,'YDir','reverse','YLim', [1 handles.numberOfRows], 'XLim', [1 handles.numberOfColumns], 'CLim',[handles.currentColorMapMinValue handles.currentColorMapMaxValue])
colormap(gray)

uistack(copyHandle(1),'bottom');


% --- Executes on button press in FradioButton.
function FradioButton_Callback(hObject, eventdata, handles)
% hObject    handle to FradioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FradioButton
set(hObject,'Value',1);
set(handles.deltaFoverFoRadioButton,'Value',0)
set(handles.ratiometricRadioButton,'Value',0)
handles.analysisMode='F';
guidata(handles.figure1,handles) ;
% --- Executes on button press in deltaFoverFoRadioButton.
function deltaFoverFoRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to deltaFoverFoRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deltaFoverFoRadioButton
set(hObject,'Value',1);
set(handles.FradioButton,'Value',0);
set(handles.ratiometricRadioButton,'Value',0)
handles.analysisMode='deltaF/F';
guidata(handles.figure1,handles) ;

%--------------------------------------------------------------------------







% --- Executes on selection change in listboxROIs.
function listboxROIs_Callback(hObject, eventdata, handles)
% hObject    handle to listboxROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listboxROIs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxROIs


% --- Executes during object creation, after setting all properties.
function listboxROIs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddNewROI.
function AddNewROI_Callback(hObject, eventdata, handles)
% hObject    handle to AddNewROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.roiNumber=handles.roiNumber+1;
set(handles.listboxROIs,'String',num2str([1:handles.roiNumber]'),...
    'Value',handles.roiNumber)
guidata(handles.figure1,handles)

% --- Executes on button press in DeleteROI.
function DeleteROI_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roiActive=get(handles.listboxROIs,'Value');
handles.roiNumber=handles.roiNumber-1;
delete(handles.roiPlotHandle{:},handles.roiTextHandle(:) );
dummy.roipolyCoord=handles.roipolyCoord;
dummy.ROIPixelList=handles.ROIPixelList;

if roiActive<=handles.roiNumber& roiActive~=1;
    i=[1:roiActive-1 roiActive+1:handles.roiNumber+1];
    set(handles.listboxROIs,'String',num2str([1:handles.roiNumber]'),...
        'Value',roiActive-1);
    
elseif roiActive==handles.roiNumber+1
    i=[1:handles.roiNumber];
    set(handles.listboxROIs,'String',num2str([1:handles.roiNumber]'),...
        'Value',roiActive-1);
else
    i=[2:handles.roiNumber+1];
    set(handles.listboxROIs,'String',num2str([1:handles.roiNumber]'),...
        'Value',roiActive);
end

handles.roipolyCoord=[];
handles.roiPlotHandle=[];
handles.roiTextHandle=[];
handles.ROIPixelList=[];


for j=1:handles.roiNumber;
    roiToUse=i(j);
    handles.roipolyCoord{j}=dummy.roipolyCoord{roiToUse};
    handles.ROIPixelList{j}=dummy.ROIPixelList{roiToUse};
    ho=plot(handles.roipolyCoord{j}.outer(:,1),handles.roipolyCoord{j}.outer(:,2),'Color',handles.colorlist(j,:),'LineWidth',2);
    if ~isempty(handles.roipolyCoord{j}.inner)
        hi=plot(handles.roipolyCoord{j}.inner(:,1),handles.roipolyCoord{j}.inner(:,2),'Color',handles.colorlist(j,:),'LineWidth',2);
    else
        hi=[];
    end
    handles.roiPlotHandle{j}=[ho hi];
    handles.roiTextHandle(j)=text(mean(handles.roipolyCoord{j}.outer(:,1)), mean(handles.roipolyCoord{j}.outer(:,2)), num2str(j),'Color',handles.colorlist(j,:),'FontWeight','bold','FontSize',8);
end
clear('dummy');
guidata(handles.figure1,handles)


% --- Executes on button press in translateROIs.
function translateROIs_Callback(hObject, eventdata, handles)
% hObject    handle to translateROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roiActive=get(handles.listboxROIs,'Value');
hAx = handles.axes1;
% CurrentPoint = get(hAx,'CurrentPoint');
hObj = handles.roiPlotHandle(roiActive);
curr_pt = get(hAx,'CurrentPoint');
% x = get(hObj,'XData');
% y = get(hObj,'YData');
%
%
% LastPoint = CurrentPoint;
% CurrentPoint = get(hAx,'CurrentPoint');
%
%
% ROIsHandles=findobj('Tag','ROI');
%
%         %Calculate Difference
%         dx = CurrentPoint(1,1) - LastPoint(1,1);
%         dy = CurrentPoint(1,2) - LastPoint(1,2);
% %         x = x + dx;
% %         y = y + dy;
%         %Update position on Plot
%      for i=1:length(ROIsHandles)
%          xy{i}(:,1)=get(ROIsHandles(i),'XData')';
%          xy{i}(:,2)=get(ROIsHandles(i),'YData')';
%         set(ROIsHandles(i),'XData',xy{i}(:,1)+dx,'YData',xy{i}(:,2)+dy);
% %         set(handles.roiTextHandle(i),'Position', [mean(xy{i}(:,1)+dx) mean(xy{i}(:,2)+dy)]);
%      end
% %   handles.roipolyCoord=xy;
%
% ud.CurrentPoint = CurrentPoint;


% set(hObj,'ButtonDownFcn','Leo_moverotateobj(''BtnDown'')')




% --- Executes on button press in lrftButton.
function lrftButton_Callback(hObject, eventdata, handles)
% hObject    handle to lrftButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i=1:length(handles.roiPlotHandle)
    outer{i}(:,1)=get(handles.roiPlotHandle{i}(1),'XData')-1';
    outer{i}(:,2)=get(handles.roiPlotHandle{i}(1),'YData')';
    
    set(handles.roiPlotHandle{i}(1),'XData',outer{i}(:,1),'YData',outer{i}(:,2));
    set(handles.roiTextHandle(i),'Position', [mean(outer{i}(:,1)) mean(outer{i}(:,2))]);
    O=find(poly2mask(outer{i}(:,1),outer{i}(:,2),handles.numberOfRows,handles.numberOfColumns));
    if length(handles.roiPlotHandle{i})>1;
        inner{i}(:,1)=get(handles.roiPlotHandle{i}(2),'XData')-1';
        inner{i}(:,2)=get(handles.roiPlotHandle{i}(2),'YData')';
        set(handles.roiPlotHandle{i}(2),'XData',inner{i}(:,1),'YData',inner{i}(:,2));
        handles.roipolyCoord{i}.inner=inner{i};
        I=find(poly2mask(inner{1}(1,:),inner{1}(2,:),handles.numberOfRows,handles.numberOfColumns));
        handles.ROIPixelList{i}.nuclear=I;
        
        handles.ROIPixelList{i}.outer=setdiff(O,I);
    else
        handles.ROIPixelList{i}.outer=O;
    end
    handles.roipolyCoord{i}.outer=outer{i};
end
guidata(handles.figure1,handles)
% --- Executes on button press in upButton.
function upButton_Callback(hObject, eventdata, handles)
% hObject    handle to upButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%ROIsHandles=findobj('Tag','ROI');
for i=1:length(handles.roiPlotHandle)
    outer{i}(:,1)=get(handles.roiPlotHandle{i}(1),'XData')';
    outer{i}(:,2)=get(handles.roiPlotHandle{i}(1),'YData')-1';
    set(handles.roiPlotHandle{i}(1),'XData',outer{i}(:,1),'YData',outer{i}(:,2));
    set(handles.roiTextHandle(i),'Position', [mean(outer{i}(:,1)) mean(outer{i}(:,2))]);
    O=find(poly2mask(outer{i}(:,1),outer{i}(:,2),handles.numberOfRows,handles.numberOfColumns));
    
    if length(handles.roiPlotHandle{i})>1;
        inner{i}(:,1)=get(handles.roiPlotHandle{i}(2),'XData')';
        inner{i}(:,2)=get(handles.roiPlotHandle{i}(2),'YData')-1';
        set(handles.roiPlotHandle{i}(2),'XData',inner{i}(:,1),'YData',inner{i}(:,2));
        handles.roipolyCoord{i}.inner=inner{i};
        I=find(poly2mask(inner{1}(1,:),inner{1}(2,:),handles.numberOfRows,handles.numberOfColumns));
        handles.ROIPixelList{i}.nuclear=I;
        handles.ROIPixelList{i}.outer=setdiff(O,I);
    else
        handles.ROIPixelList{i}.outer=O;
    end
    
    handles.roipolyCoord{i}.outer=outer{i};
    
    
end
guidata(handles.figure1,handles)
% --- Executes on button press in rightButton.
function rightButton_Callback(hObject, eventdata, handles)
% hObject    handle to rightButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i=1:length(handles.roiPlotHandle)
    outer{i}(:,1)=get(handles.roiPlotHandle{i}(1),'XData')+1';
    outer{i}(:,2)=get(handles.roiPlotHandle{i}(1),'YData')';
    set(handles.roiPlotHandle{i}(1),'XData',outer{i}(:,1),'YData',outer{i}(:,2));
    set(handles.roiTextHandle(i),'Position', [mean(outer{i}(:,1)) mean(outer{i}(:,2))]);
    O=find(poly2mask(outer{i}(:,1),outer{i}(:,2),handles.numberOfRows,handles.numberOfColumns));
    
    if length(handles.roiPlotHandle{i})>1;
        inner{i}(:,1)=get(handles.roiPlotHandle{i}(2),'XData')+1';
        inner{i}(:,2)=get(handles.roiPlotHandle{i}(2),'YData')';
        set(handles.roiPlotHandle{i}(2),'XData',inner{i}(:,1),'YData',inner{i}(:,2));
        handles.roipolyCoord{i}.inner=inner{i};
        I=find(poly2mask(inner{1}(1,:),inner{1}(2,:),handles.numberOfRows,handles.numberOfColumns));
        handles.ROIPixelList{i}.nuclear=I;
        handles.ROIPixelList{i}.outer=setdiff(O,I);
    else
        handles.ROIPixelList{i}.outer=O;
    end

    handles.roipolyCoord{i}.outer=outer{i};
    
    
end
guidata(handles.figure1,handles)
% --- Executes on button press in downButton.
function downButton_Callback(hObject, eventdata, handles)
% hObject    handle to downButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i=1:length(handles.roiPlotHandle)
    outer{i}(:,1)=get(handles.roiPlotHandle{i}(1),'XData')';
    outer{i}(:,2)=get(handles.roiPlotHandle{i}(1),'YData')+1';
    set(handles.roiPlotHandle{i}(1),'XData',outer{i}(:,1),'YData',outer{i}(:,2));
    set(handles.roiTextHandle(i),'Position', [mean(outer{i}(:,1)) mean(outer{i}(:,2))]);
    O=find(poly2mask(outer{i}(:,1),outer{i}(:,2),handles.numberOfRows,handles.numberOfColumns));
    if length(handles.roiPlotHandle{i})>1;
        inner{i}(:,1)=get(handles.roiPlotHandle{i}(2),'XData')';
        inner{i}(:,2)=get(handles.roiPlotHandle{i}(2),'YData')+1';
        set(handles.roiPlotHandle{i}(2),'XData',inner{i}(:,1),'YData',inner{i}(:,2));
        handles.roipolyCoord{i}.inner=inner{i};
        I=find(poly2mask(inner{1}(1,:),inner{1}(2,:),handles.numberOfRows,handles.numberOfColumns));
        handles.ROIPixelList{i}.nuclear=I;
        handles.ROIPixelList{i}.outer=setdiff(O,I);
    else
        handles.ROIPixelList{i}.outer=O;
    end
    handles.roipolyCoord{i}.outer=outer{i};
end
guidata(handles.figure1,handles)


% --- Executes on button press in ratiometricRadioButton.
function ratiometricRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to ratiometricRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ratiometricRadioButton

set(hObject,'Value',1);
set(handles.FradioButton,'Value',0);
set(handles.deltaFoverFoRadioButton,'Value',0)
handles.analysisMode='deltaF/F/R';
guidata(handles.figure1,handles) ;


% --- Executes on button press in displayGreenRadioButton.
function displayGreenRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to displayGreenRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayGreenRadioButton
set(hObject,'Value',1);
set(handles.displayRedRadioButton,'Value',0);
handles.displayChannel='G';
guidata(handles.figure1,handles) ;

% --- Executes on button press in displayRedRadioButton.
function displayRedRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to displayRedRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayRedRadioButton

set(hObject,'Value',1);
set(handles.displayGreenRadioButton,'Value',0);
handles.displayChannel='R';
guidata(handles.figure1,handles) ;



function startingTrialBox_Callback(hObject, eventdata, handles)
% hObject    handle to startingTrialBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startingTrialBox as text
%        str2double(get(hObject,'String')) returns contents of startingTrialBox as a double


% --- Executes during object creation, after setting all properties.
function startingTrialBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startingTrialBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endingTrialBox_Callback(hObject, eventdata, handles)
% hObject    handle to endingTrialBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endingTrialBox as text
%        str2double(get(hObject,'String')) returns contents of endingTrialBox as a double


% --- Executes during object creation, after setting all properties.
function endingTrialBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endingTrialBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in analyzeAcrossTrialsButton.
function analyzeAcrossTrialsButton_Callback(hObject, eventdata, handles)
% hObject    handle to analyzeAcrossTrialsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'channels')
    handles.channels=[];
end
startString=get(handles.startingTrialBox,'String');
endingString=get(handles.endingTrialBox,'String');
startTrialNum=str2num(startString);
endingTrialNum=str2num(endingString);
%handles.lastF=zeros(handles.numberOfFrames);
%handles.FAcrossTrials=zeros(handles.roiNumber, handles.numberOfFrames, endingTrialNum-startTrialNum+1);
clear('handles.trials');
tic

% % get header from file from ScanImage to compensate for lack of header info
% % in EP files
% [f2,p2]=uigetfile('*.tif','Choose file with intact header');
% 
% % insert all the necessary header info here
% [header,~,~,~] = scanimage.util.opentif(fullfile(p2,f2));
H.scanZoomFactor = handles.state.SI.hRoiManager.scanZoomFactor;
H.scanLinePeriod = handles.state.SI.hRoiManager.linePeriod;
H.scanFrameRate = handles.state.SI.hRoiManager.scanFrameRate;
H.beamPowers = handles.state.SI.hBeams.powers;
H.fastZEnable = handles.state.SI.hFastZ.enable;
H.stackZStepSize = handles.state.SI.hStackManager.stackZStepSize;
H.stackNumSlices = 1;
H.channelsActive = handles.state.SI.hChannels.channelsActive;
H.date = handles.state.date;
handles.timePerLine=H.scanLinePeriod;

handles.state=H;
numPlanes=handles.state.stackNumSlices;
handles.numberOfPlanes=numPlanes;

handles.startTrialNum=startTrialNum;
handles.endingTrialNum=endingTrialNum;

%loop through trials
[F,P]=uigetfile('*.tif','Pick first trial TIFF');
firstTrial_SI=split(F,'_'); firstTrial_SI=str2double(firstTrial_SI{end}(1:end-4)); 
trialIndexes=firstTrial_SI + (0:endingTrialNum-1);
for i=startTrialNum:endingTrialNum;
    f=[F(1:end-7) sprintf('%03d',trialIndexes(i)) '.tif'];
    handles.currentFile=f;
    handles.currentTrialNumber=i;
    
    

    if exist(f,'file')
        
        [~,Data]=opentif_SI5a(f);
        disp([f ' loaded'])
        Data(Data==0)=NaN;
        handles.channels.G = Data;
        
        
        [r c z]=size(handles.channels.G);
        handles.numberOfFrames=z;
        if get(handles.loadAvgFilterCheck,'Value')==1
            filter=fspecial('average',[2 2]);
            disp('filtering....')
            fImage=imfilter(handles.channels.G, filter, 'conv');
            handles.channels.G=fImage;
        end
        
        
        guidata(handles.figure1,handles)
        AnalyzeButton_Callback(hObject, eventdata, handles)
        handles=guidata(handles.figure1);
        handles.FAcrossTrials{i}=handles.lastF;
        handles.F_Neuropil_AcrossTrials{i}=handles.lastFNeuropil;% rows= ROI , column= frame, z=trial    handles.trials(i,:)=filename;
        handles.backgroundTraces{i,:}=handles.lastBackgroundTrace;
    end
end
toc

guidata(handles.figure1,handles) ;
% stuff to include in the header of the Analysis results:
FAcrossTrials.data=handles.FAcrossTrials;
FAcrossTrials.NeuropilData=handles.F_Neuropil_AcrossTrials;
FAcrossTrials.header.directory=handles.currentDirectory;
FAcrossTrials.header.filename=handles.currentFile;
FAcrossTrials.header.startTrialNumber=startTrialNum;
FAcrossTrials.header.endTrialNumber=endingTrialNum;
FAcrossTrials.header.timePerLine=handles.timePerLine;
FAcrossTrials.header.numberOfRows=handles.numberOfRows;
FAcrossTrials.header.analysisMode=handles.analysisMode;
FAcrossTrials.header.roipolyCoord=handles.roipolyCoord;
FAcrossTrials.header.backgroundpolyCoord=handles.backgroundpolyCoord;
FAcrossTrials.header.background=handles.background;
FAcrossTrials.header.backgroundTraces=handles.backgroundTraces;

eval([['FAcrossTrials_' handles.currentFile(1:end-10)],'=FAcrossTrials']);
save(['FAcrossTrials_' handles.currentFile(1:end-10) '.mat'],['FAcrossTrials_' handles.currentFile(1:end-10)])
assignin('base', ['FAcrossTrials_' handles.currentFile(1:end-10)], FAcrossTrials);




function ROItoPlotValue_Callback(hObject, eventdata, handles)
% hObject    handle to ROItoPlotValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ROItoPlotValue as text
%        str2double(get(hObject,'String')) returns contents of ROItoPlotValue as a double


% --- Executes during object creation, after setting all properties.
function ROItoPlotValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROItoPlotValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotROIDataButton.
function plotROIDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to plotROIDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ROItoPlot=str2num(get(handles.ROItoPlotValue,'String'));
figure(111);
timePerLine = handles.timePerLine; % in s
deltaT = handles.numberOfRows * timePerLine*handles.numberOfPlanes; % in ms, time per frame
yFactor=20;    %<-----Offset between traces for the plot

Data=handles.FAcrossTrials;
[r c]=size(Data);
offsetVector=yFactor*(0:c-1)';

sizeToPad=65;
%sizeToPad=min(cellfun(@length,Data));
offsetArray=repmat(offsetVector, [1 sizeToPad ]);
for i=1:length(Data);
    
    ROIDataPadded(i,:)=Data{i}(ROItoPlot,1:sizeToPad);
end
ROIplot=ROIDataPadded+offsetArray;

plot(deltaT:deltaT:sizeToPad*deltaT,ROIplot', 'Color', handles.colorlist(ROItoPlot,:))

hax=gca;
set(get(hax,'XLabel'),'String', 'Time (s)');
set(get(hax,'YLabel'),'String', [handles.analysisMode '--Trial #']);
set(hax,'XLim',[0 (handles.numberOfFrames/handles.numberOfPlanes)*deltaT]);
set(hax,'YLim',[-0.5 c+600]);
title(['ROI ' num2str(ROItoPlot)]);
figure(112)
plot(deltaT:deltaT:sizeToPad*deltaT,mean(ROIDataPadded));

hax=gca;
set(get(hax,'XLabel'),'String', 'Time (s)');
set(get(hax,'YLabel'),'String', [handles.analysisMode '--Trial #']);
% set(hax,'XLim',[0 (handles.numberOfFrames/handles.numberOfPlanes)*deltaT]);
% set(hax,'YLim',[-0.5 c+600]);
figure(122)
plot(reshape(ROIDataPadded',[1 c*sizeToPad]))


% --- Executes on button press in saveROIbutton.
function saveROIbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saveROIbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
savedROIs.backgroundpolyCoord=handles.backgroundpolyCoord;
savedROIs.roipolyCoord=handles.roipolyCoord;
savedROIs.ROIPixelList=handles.ROIPixelList;

uisave('savedROIs')
%save(['savedROIs_' handles.currentDirectory(end-17:end-1)],'savedROIs');
disp(['saving.. savedROIs ' handles.currentDirectory])




% --- Executes on button press in loadROIsButton.
function loadROIsButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadROIsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f p]=uigetfile('Pick saved ROIs');
loadedROIs=load([p f]);
try;
    delete(handles.roiPlotHandle(:),handles.roiTextHandle(:), handles.backgroundHandle )
catch
    handles.roiPlotHandle=[];
    handles.roiTextHandle=[];
end
handles.roiPlotHandle=[];
handles.roiTextHandle=[];
handles.roipolyCoord=[];
handles.backgroundpolyCoord=[];
handles.ROIPixelList=[];
[r c]=size(loadedROIs.savedROIs.roipolyCoord);

handles.roiNumber=c;
set(handles.listboxROIs,'String',num2str([1:handles.roiNumber]'),...
    'Value',1);

handles.roipolyCoord=loadedROIs.savedROIs.roipolyCoord;
handles.backgroundpolyCoord=loadedROIs.savedROIs.backgroundpolyCoord;
handles.ROIPixelList=loadedROIs.savedROIs.ROIPixelList;
%calculate background;


%plot:
hold on
handles.backgroundHandle=plot(handles.backgroundpolyCoord(:,1), handles.backgroundpolyCoord(:,2),'w-');


for i=1:handles.roiNumber;
    
    
    ho=plot(handles.roipolyCoord{i}.outer(:,1),handles.roipolyCoord{i}.outer(:,2),'Color',handles.colorlist(i,:),'LineWidth',2);
    if ~isempty(handles.roipolyCoord{i}.inner)
        hi=plot(handles.roipolyCoord{i}.inner(:,1),handles.roipolyCoord{i}.inner(:,2),'Color',handles.colorlist(i,:),'LineWidth',2);
    else
        hi=[];
    end
    
    handles.roiPlotHandle{i}=[ho hi];
    handles.roiTextHandle(i)=text(mean(handles.roipolyCoord{i}.outer(:,1)), mean(handles.roipolyCoord{i}.outer(:,2)), num2str(i),'Color',handles.colorlist(i,:),'FontWeight','bold','FontSize',8);
end

guidata(handles.figure1,handles) ;


% --- Executes on button press in loadFowardButton.
function loadFowardButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadFowardButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadFowardButton
filename=handles.currentFile;
number=str2num(filename(end-6:end-4));
filename=[handles.currentFile(1:end-7) sprintf('%03d',number+1) '.tif'];
f=filename;
handles.currentFile=filename;
if isfield(handles,'imageHandle')
    delete(handles.imageHandle)
    
end
if isfield(handles,'ccimage');
    handles=rmfield(handles, 'ccimage');
end

%load file
[H,Data] = load_SI2019a_EP(handles.currentDirectory,f);
handles.timePerLine=H.scanLinePeriod;
Data(Data==0)=NaN;
handles.state=H;
handles.channels.G=squeeze(Data(:,:,1,:));

if get(handles.loadAvgFilterCheck,'Value')==1
    filter=fspecial('average',[2 2]);
    disp('filtering....')
    fImage=imfilter(handles.channels.G, filter, 'conv');
    handles.channels.G=fImage;
end
s=size(Data);
if s(3)==2
    handles.channels.R=squeeze(Data(:,:,2,:));
    if get(handles.loadAvgFilterCheck,'Value')==1
        filter=fspecial('average',[2 2]);
        disp('filtering....')
        fImage=imfilter(handles.channels.R, filter, 'conv');
        handles.channels.R=fImage;
    end
end

[r c z]=size(handles.channels.G);
handles.numberOfFrames=z;
handles.numberOfRows=r;
handles.numberOfColumns=c;
set(handles.figure1,'Name',handles.currentFile);
set(handles.slider1,'Min',1,'Max',z,'Value',1,'SliderStep',[ 1/z 10/z]);
set(handles.FrameCounter,'String',[num2str(1) '/' num2str(z)]);
set(handles.zoomFactorBox,'String', ['Zoom: ' num2str(handles.state.scanZoomFactor)]);
set(handles.DateBox,'String', ['Date: '  handles.state.date]);
set(handles.PowerBox,'String', ['Power: ' num2str(handles.state.beamPowers)]);
handles.timePerLine=handles.state.scanLinePeriod; %in sec
guidata(handles.figure1,handles) ;


MaxProjection_Callback(hObject, eventdata, handles);



% --- Executes on button press in loadBackwardsButton.
function loadBackwardsButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadBackwardsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadBackwardsButton
filename=handles.currentFile;
number=str2num(filename(end-6:end-4));
filename=[handles.currentFile(1:end-7) sprintf('%03d',number-1) '.tif'];
f=filename;
handles.currentFile=filename;
if isfield(handles,'imageHandle')
    delete(handles.imageHandle)
    
end
if isfield(handles,'ccimage');
    handles=rmfield(handles, 'ccimage');
end


[H,Data] = load_SI2019a_EP(handles.currentDirectory,f);

Data(Data==0)=NaN;
handles.state=H;
handles.channels.G=squeeze(Data(:,:,1,:));

if get(handles.loadAvgFilterCheck,'Value')==1
    filter=fspecial('average',[2 2]);
    disp('filtering....')
    fImage=imfilter(handles.channels.G, filter, 'conv');
    handles.channels.G=fImage;
end
s=size(Data);
if s(3)==2
    handles.channels.R=squeeze(Data(:,:,2,:));
    if get(handles.loadAvgFilterCheck,'Value')==1
        filter=fspecial('average',[2 2]);
        disp('filtering....')
        fImage=imfilter(handles.channels.R, filter, 'conv');
        handles.channels.R=fImage;
    end
end

[r c z]=size(handles.channels.G);
handles.numberOfFrames=z;
handles.numberOfRows=r;
handles.numberOfColumns=c;
set(handles.figure1,'Name',handles.currentFile);
set(handles.slider1,'Min',1,'Max',z,'Value',1,'SliderStep',[ 1/z 10/z]);
set(handles.FrameCounter,'String',[num2str(1) '/' num2str(z)]);
set(handles.zoomFactorBox,'String', ['Zoom: ' num2str(handles.state.scanZoomFactor)]);
set(handles.DateBox,'String', ['Date: '  handles.state.date]);
set(handles.PowerBox,'String', ['Power: ' num2str(handles.state.beamPowers)]);
handles.timePerLine=handles.state.scanLinePeriod; %in sec
guidata(handles.figure1,handles) ;


MaxProjection_Callback(hObject, eventdata, handles);




% --- Executes on button press in ROItoPlotBackButton.
function ROItoPlotBackButton_Callback(hObject, eventdata, handles)
% hObject    handle to ROItoPlotBackButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ROItoPlotBackButton
ROItoPlot=str2num(get(handles.ROItoPlotValue,'String'));
ROItoPlot=ROItoPlot-1;
set(handles.ROItoPlotValue,'String',num2str(ROItoPlot));
guidata(handles.figure1,handles) ;
plotROIDataButton_Callback(hObject, eventdata, handles)
%plotTrialButton_Callback(hObject, eventdata, handles)

% --- Executes on button press in ROItoPlotForwardButton.
function ROItoPlotForwardButton_Callback(hObject, eventdata, handles)
% hObject    handle to ROItoPlotForwardButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ROItoPlotForwardButton

ROItoPlot=str2num(get(handles.ROItoPlotValue,'String'));
ROItoPlot=ROItoPlot+1;
set(handles.ROItoPlotValue,'String',num2str(ROItoPlot));
guidata(handles.figure1,handles) ;
plotROIDataButton_Callback(hObject, eventdata, handles)
%plotTrialButton_Callback(hObject, eventdata, handles)



% --- Executes on button press in loadAvgFilterCheck.
function loadAvgFilterCheck_Callback(hObject, eventdata, handles)
% hObject    handle to loadAvgFilterCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadAvgFilterCheck






% --- Executes on button press in plotTrialButton.
function plotTrialButton_Callback(hObject, eventdata, handles)
% hObject    handle to plotTrialButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

trialToPlot=str2num(get(handles.ROItoPlotValue,'String'));
figure(112);
timePerLine = handles.timePerLine; % in s
deltaT = handles.numberOfRows * timePerLine; % in ms, time per frame
yFactor=1;    %<-----Offset between traces for the plot
%[a b c]=size(handles.FAcrossTrials);
trialData=handles.FAcrossTrials{trialToPlot};
offsetVector=yFactor*(0:handles.roiNumber-1)';
offsetArray=repmat(offsetVector, [1 handles.numberOfFrames/handles.numberOfPlanes ]);
trialPlot=trialData+offsetArray;
set(gcf,'DefaultAxesColorOrder',handles.colorlist);
plot(deltaT:deltaT:(handles.numberOfFrames/handles.numberOfPlanes)*deltaT,trialPlot')

hax=gca;
title(['Trial# ' num2str(trialToPlot)]);
set(get(hax,'XLabel'),'String', 'Time (s)');
set(get(hax,'YLabel'),'String', [handles.analysisMode '-- ROI #']);
set(hax,'XLim',[0 handles.numberOfFrames*deltaT]);
set(hax,'YLim',[-0.5 1*handles.roiNumber+2]);



% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)


%
% try
%     delete(handles.imageHandle)
% catch
% end
%
% switch handles.displayChannel
%     case 'G'
%         activeChannel=handles.channels.G;
%     case 'R'
%         activeChannel=handles.channels.R;
% end
%
% currentFrame = round(get(handles.slider1,'Value'));
%
% frame=currentFrame+ eventdata.VerticalScrollCount;
%
% if frame<=0
%     frame=1;
% end
%
% if frame>handles.numberOfFrames;
%     frame=handles.numberOfFrames;
% end
%
% im=imagesc(activeChannel(:,:,frame),[handles.currentColorMapMinValue handles.currentColorMapMaxValue]);
% handles.imageHandle=im;
% uistack(handles.imageHandle,'bottom')
% set(handles.slider1,'Value', frame);
% set(handles.FrameCounter,'String',[num2str(frame) '/' num2str(handles.numberOfFrames) ' | ' num2str(handles.numberOfRows * handles.timePerLine*frame) 's']);
% guidata(handles.figure1,handles) ;




% --- Executes on selection change in targetGenerationMethod.
function targetGenerationMethod_Callback(hObject, eventdata, handles)
% hObject    handle to targetGenerationMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns targetGenerationMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from targetGenerationMethod

value=get(hObject,'Value');
string=get(hObject,'String');
handle.targetGenerationMethodString=string(value,:);
guidata(handles.figure1,handles) ;

% --- Executes during object creation, after setting all properties.
function targetGenerationMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetGenerationMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in generateTarget.
function generateTarget_Callback(hObject, eventdata, handles)
% hObject    handle to generateTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value=get(handles.targetGenerationMethod,'Value');
string=get(handles.targetGenerationMethod,'String');
handles.targetGenerationMethodString=string(value,:);
saveName=get(handles.targetSaveName,'String');
currentFrame = round(get(handles.slider1,'Value'));

[pathstr, name, ~]= fileparts([handles.currentDirectory, handles.currentFile]);

switch handles.displayChannel
    case 'G'
        channel=1;
        if handles.numberOfPlanes>1
            activeChannel=handles.channels.G{handles.planeOfInterest};
        else
            activeChannel=handles.channels.G;
        end
    case 'R'
        channel=2;
        if handles.numberOfPlanes>1
            activeChannel=handles.channels.R{handles.planeOfInterest};
        else
            activeChannel=handles.channels.R;
        end
end


switch handles.targetGenerationMethodString{1};
    case 'Trial Reg on Frame'
        
        
        target=activeChannel(:,:,currentFrame);
        %for red and green channel registration
%         Leo_image_registration_MultPeak_singleTrial_padding_singlePlane(target,fullfile(handles.currentDirectory, handles.currentFile), [saveName 'on_frame_' num2str(currentFrame) '_'],[1 2],handles.planeOfInterest,handles.numberOfPlanes);
        
        %for only green channel registration
        Leo_image_registration_MultPeak_singleTrial_padding_singlePlane(target,fullfile(handles.currentDirectory, handles.currentFile), [saveName 'on_frame_' num2str(currentFrame) '_'],channel,handles.planeOfInterest,handles.numberOfPlanes);
        
%         Marina_image_registration_single_trial_padding_singlePlane(target,fullfile(handles.currentDirectory, handles.currentFile), [saveName 'on_frame_' num2str(currentFrame) '_'],channel,handles.planeOfInterest,handles.numberOfPlanes,1);
        %Leo_image_registration_single_trial_single_plane(target,fullfile(handles.currentDirectory, handles.currentFile), [saveName 'on_frame_' num2str(currentFrame) '_'],channel,handles.planeOfInterest,handles.numberOfPlanes);
%         Leo_image_registration_single_trial_padding_singlePlane(target,fullfile(handles.currentDirectory, handles.currentFile), [saveName 'on_frame_' num2str(currentFrame) '_'],channel,handles.planeOfInterest,handles.numberOfPlanes);
        
        %         [regTargetTrial  header]=loadScanImageMovieSingleChannel(fullfile(handles.currentDirectory, [saveName 'on_frame_' num2str(currentFrame) '_' num2str(channel) name '.tif' ]),pathstr,1);
        %         maxProjOfRegisteredTrial=mean(im_mov_avg(regTargetTrial,5),3);
        %         imwrite(uint8(maxProjOfRegisteredTrial), fullfile(handles.currentDirectory, [saveName num2str(currentFrame) 'plane' num2str(handles.planeOfInterest) '_mean_movAvg.tif']), 'tif');
        
    case 'Trial Max.'
        
        
%         regTargetTrial=im_mov_avg(activeChannel,5);
%         maxProjOfRegisteredTrial=mean(im_mov_avg(regTargetTrial,5),3);
        
        %-- gabi april 2021
        maxProjOfRegisteredTrial=nanmean(activeChannel,3);
        maxProjOfRegisteredTrial=maxProjOfRegisteredTrial-(min(maxProjOfRegisteredTrial(:))*1.1);
        %--

        imwrite(uint16(maxProjOfRegisteredTrial), fullfile(handles.currentDirectory, [saveName num2str(currentFrame) '_mean_movAvg.tif']), 'tif');
    case 'Frame'
        
        
        target=activeChannel(:,:,currentFrame);
        
        
        
        imwrite(uint8(target), fullfile(handles.currentDirectory, [saveName num2str(currentFrame) '.tif']), 'tif');
end

function targetSaveName_Callback(hObject, eventdata, handles)
% hObject    handle to targetSaveName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetSaveName as text
%        str2double(get(hObject,'String')) returns contents of targetSaveName as a double


% --- Executes during object creation, after setting all properties.
function targetSaveName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetSaveName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function currentColorMapMin_Callback(hObject, eventdata, handles)
% hObject    handle to currentColorMapMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentColorMapMin as text
%        str2double(get(hObject,'String')) returns contents of currentColorMapMin as a double
handles.currentColorMapMaxValue=str2double(get(handles.currentColorMapMax,'String'));
handles.currentColorMapMinValue=str2double(get(hObject,'String'));
set(handles.ColormapMin, 'Value', handles.currentColorMapMinValue);
guidata(handles.figure1,handles) ;
set(handles.axes1,'CLim', [handles.currentColorMapMinValue handles.currentColorMapMaxValue]);


% --- Executes during object creation, after setting all properties.
function currentColorMapMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentColorMapMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function currentColorMapMax_Callback(hObject, eventdata, handles)
% hObject    handle to currentColorMapMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentColorMapMax as text
%        str2double(get(hObject,'String')) returns contents of
%        currentColorMapMax as a double
handles.currentColorMapMaxValue=str2double(get(hObject,'String'));
handles.currentColorMapMinValue=str2double(get(handles.currentColorMapMin,'String'));
set(handles.ColormapMax, 'Value', handles.currentColorMapMaxValue);
guidata(handles.figure1,handles) ;
set(handles.axes1,'CLim', [handles.currentColorMapMinValue handles.currentColorMapMaxValue]);

% --- Executes during object creation, after setting all properties.
function currentColorMapMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentColorMapMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in FindROIButton.
function FindROIButton_Callback(hObject, eventdata, handles)
% hObject    handle to FindROIButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roiActive=get(handles.listboxROIs,'Value');
currentLineWidth=get(handles.roiPlotHandle{roiActive}(1), 'LineWidth');
if currentLineWidth==2
    set(handles.roiPlotHandle{roiActive}(1), 'LineWidth',4)
else
    set(handles.roiPlotHandle{roiActive}(1), 'LineWidth',2)
end








% --- Executes on button press in SDProjection.
function SDProjection_Callback(hObject, eventdata, handles)
% hObject    handle to SDProjection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'Name',handles.currentFile);
set(handles.ColormapMax,'Min',0, 'Max',max(max(max(handles.channels.G,[],3))), 'Value', max(max(max(handles.channels.G,[],3))),...
    'SliderStep', [ 1/double(max(handles.channels.G(:))) 10/double(max(handles.channels.G(:)))])%
set(handles.ColormapMin,'Min',0, 'Max',max(max(max(handles.channels.G,[],3))), 'Value', 0,...
    'SliderStep', [ 1/double(max(handles.channels.G(:))) 10/double(max(handles.channels.G(:)))])%
set(handles.MaxValueText, 'String', [ '/' num2str(max(max(max(handles.channels.G,[],3))))]);
set(handles.MinValueText, 'String', ['/' num2str(min(min(min(handles.channels.G,[],3))))]);
try
    delete(handles.imageHandle)
    
catch
end

switch handles.displayChannel
    case 'G'
        SDim=nanstd(double(handles.channels.G),[],3);%
    case 'R'
        SDim=nanstd(double(handles.channels.R),[],3);%
    case 'ProjectAcrossTrials'
        SDim=nanstd(double(handles.trialsProjection),[],3);
        
end
axes(handles.axes1);
colormapMax=str2double(get(handles.currentColorMapMax,'String'));
colormapMin=str2double(get(handles.currentColorMapMin,'String'));

im = showImage(SDim,handles);
colormap(gray);
handles.imageHandle=im;
uistack(handles.imageHandle,'bottom')

guidata(handles.figure1,handles) ;


% --- Executes on button press in AvgButton.
function AvgButton_Callback(hObject, eventdata, handles)
% hObject    handle to AvgButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'Name',handles.currentFile);
if handles.numberOfPlanes == 1
    set(handles.ColormapMax,'Min',0, 'Max',max(max(max(handles.channels.G,[],3))), 'Value', 0,...
        'SliderStep', [ 1/double(max(handles.channels.G(:))) 10/double(max(handles.channels.G(:)))])%
    set(handles.ColormapMin,'Min',min(min(min(handles.channels.G,[],3))), 'Max',max(max(max(handles.channels.G,[],3))), 'Value',200,...
        'SliderStep', [ 1/double(max(handles.channels.G(:))) 10/double(max(handles.channels.G(:)))])%
    set(handles.MaxValueText, 'String', [ '/' num2str(max(max(max(handles.channels.G,[],3))))]);
    set(handles.MinValueText, 'String', ['/' num2str(min(min(min(handles.channels.G,[],3))))]);
else
    set(handles.ColormapMax,'Min',0, 'Max',max(max(max(handles.channels.G{handles.planeOfInterest},[],3))), 'Value', 0,...
        'SliderStep', [ 1/double(max(handles.channels.G{handles.planeOfInterest}(:))) 10/double(max(handles.channels.G{handles.planeOfInterest}(:)))])%
    set(handles.ColormapMin,'Min',min(min(min(handles.channels.G{handles.planeOfInterest},[],3))), 'Max',max(max(max(handles.channels.G{handles.planeOfInterest},[],3))), 'Value',200,...
        'SliderStep', [ 1/double(max(handles.channels.G{handles.planeOfInterest}(:))) 10/double(max(handles.channels.G{handles.planeOfInterest}(:)))])%
    set(handles.MaxValueText, 'String', [ '/' num2str(max(max(max(handles.channels.G{handles.planeOfInterest},[],3))))]);
    set(handles.MinValueText, 'String', ['/' num2str(min(min(min(handles.channels.G{handles.planeOfInterest},[],3))))]);
end

try
    delete(handles.imageHandle)   
catch
end

switch handles.displayChannel
    case 'G'
        if handles.numberOfPlanes>1
            Avgim=nanmean(handles.channels.G{handles.planeOfInterest},3);
        else
            Avgim=nanmean(handles.channels.G,3);%
        end
    case 'R'
        
        if handles.numberOfPlanes>1
            Avgim=nanmean(handles.channels.R{handles.planeOfInterest},3);
        else
            Avgim=nanmean(handles.channels.R,3);%
        end
        %
    case 'ProjectAcrossTrials'
        Avgim=nanmean(handles.trialsProjection,3);        
end
axes(handles.axes1);
colormapMax=str2double(get(handles.currentColorMapMax,'String'));
colormapMin=str2double(get(handles.currentColorMapMin,'String'));
im = showImage(Avgim,handles);
colormap(gray);

handles.imageHandle=im;
uistack(handles.imageHandle,'bottom')

guidata(handles.figure1,handles) ;



function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double


% --- Executes during object creation, after setting all properties.
function threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AutoROIs.
function AutoROIs_Callback(hObject, eventdata, handles)
% hObject    handle to AutoROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try;
    delete(handles.roiPlotHandle{:},handles.roiTextHandle(:))
catch
    handles.roiPlotHandle=[];
    handles.roiTextHandle=[];
end
handles.roiPlotHandle=[];
handles.roiTextHandle=[];
handles.roipolyCoord=[];

currentImage=handles.imageHandle;
threshold=str2double(get(handles.threshold,'String'));
ROIminSize=str2double(get(handles.ROIminSize,'String'));
contourWidth=str2double(get(handles.contourWidth,'String'));
trimRows=str2double(get(handles.TrimRows,'String'));
trimColumns=str2double(get(handles.TrimColumns,'String'));

currentImage=get(currentImage,'CData');
currentImage(1:trimRows,:)=0;
currentImage(handles.numberOfRows-trimRows:end,:)=0;
currentImage(:,1:trimColumns)=0;
currentImage(:,handles.numberOfColumns-trimColumns:end)=0;


BW=currentImage>threshold;
BWA=bwareaopen(BW,ROIminSize);

se90 = strel('line', contourWidth, 90);
se0 = strel('line', contourWidth, 0);


%BWAF=(imdilate(imclearborder(imfill(BWA,'holes')),[ se90 se0]));
BWAF=imclearborder(imfill(BWA,'holes'));

cc=bwconncomp(BWAF,8);
handles.roiNumber=cc.NumObjects;
hold on
for i=1:cc.NumObjects;
    I0=zeros(cc.ImageSize);
    I0(cc.PixelIdxList{i})=1;
    I0dil=imdilate(I0,[ se90 se0] );
    ROIcontour=contourc(I0dil,1);
    handles.roipolyCoord{i}.outer=(ROIcontour(:,2:end))';
    handles.roipolyCoord{i}.inner=[];
    handles.ROIPixelList{i}.outer=find(I0dil);
    handles.ROIPixelList{i}.inner=[];
    
    handles.roiPlotHandle{i}=plot(handles.roipolyCoord{i}.outer(:,1),handles.roipolyCoord{i}.outer(:,2),'Color',handles.colorlist(i,:),'LineWidth',2);
    handles.roiTextHandle(i)=text(mean(handles.roipolyCoord{i}.outer(:,1))-2*length(num2str(i)), mean(handles.roipolyCoord{i}.outer(:,2)), num2str(i),'Color',handles.colorlist(i,:),'FontWeight','bold');
end



set(handles.listboxROIs,'String',num2str([1:handles.roiNumber]'),...
    'Value',1);
uistack(handles.imageHandle,'bottom')
guidata(handles.figure1,handles) ;

uistack(handles.imageHandle,'bottom')
guidata(handles.figure1,handles) ;







function contourWidth_Callback(hObject, eventdata, handles)
% hObject    handle to contourWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contourWidth as text
%        str2double(get(hObject,'String')) returns contents of contourWidth as a double


% --- Executes during object creation, after setting all properties.
function contourWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contourWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ROIminSize_Callback(hObject, eventdata, handles)
% hObject    handle to ROIminSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ROIminSize as text
%        str2double(get(hObject,'String')) returns contents of ROIminSize as a double


% --- Executes during object creation, after setting all properties.
function ROIminSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROIminSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ProjectAcrossTrials.
function ProjectAcrossTrials_Callback(hObject, eventdata, handles)
% hObject    handle to ProjectAcrossTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

startString=get(handles.startingTrialBox,'String');
endingString=get(handles.endingTrialBox,'String');
startTrialNum=str2num(startString);
endingTrialNum=str2num(endingString);
trialsProjection=zeros(handles.numberOfRows, handles.numberOfColumns, endingTrialNum-startTrialNum+1);
clear('handles.trials');
string=get(handles.ProjectionItem,'String');
projectionItem=string(get(handles.ProjectionItem,'Value'),:);

string2=get(handles.ProjecAcrossTrialsSecondStep,'String');
projectionSecondStep=string2(get(handles.ProjecAcrossTrialsSecondStep,'Value'),:);
filename=[handles.currentFile(1:end-7) sprintf('%03d',i) '.tif'];
parfor trial=startTrialNum:endingTrialNum;
    %filename=[handles.currentFile(1:end-7) sprintf('%03d',i) '.tif'];
    f=[filename(1:end-7) sprintf('%03d',trial) '.tif'];
    %handles.currentFile=filename;
    
    disp(trial)
    if exist(f,'file')
        
        [H Data]=scim_openTif(f);
        Data(Data==0)=NaN;
        
        G=squeeze(Data(:,:,1,:));
        
        
        switch projectionItem{1}
            case 'SD'
                trialsProjection(:,:,trial)=nanstd(double(G),[],3);%
            case 'Avg Proj'
                trialsProjection(:,:,trial)=nanmean(G,3);%
            case 'Mov Avg Max Proj. '
                trialsProjection(:,:,trial)=max(im_mov_avg(G,5),[],3);%
                
        end
    end
end
try
    delete(handles.imageHandle)
    
catch
end
handles.trialsProjection=trialsProjection;


axes(handles.axes1);
colormapMax=str2double(get(handles.currentColorMapMax,'String'));
colormapMin=str2double(get(handles.currentColorMapMin,'String'));
switch projectionSecondStep{1}
    case 'Max'
        im=showImage(max(handles.trialsProjection,[],3),handles);colormap(gray);
    case 'Avg'
        im=showImage(nanmean(handles.trialsProjection,3),handles);colormap(gray);
end
handles.imageHandle=im;
uistack(handles.imageHandle,'bottom')
%set sliderTrials
trial=round(get(handles.sliderTrials,'Value'));
set(handles.sliderTrials,'Min',1,'Max',endingTrialNum,'Value',1,'SliderStep',[ 1/endingTrialNum 10/endingTrialNum]);
set(handles.TrialCounter,'String',[num2str(trial) '/' endingString ]);
handles.currentFile=filename;
guidata(handles.figure1,handles) ;



assignin('base', ['ProjectionAcrossTrials_' handles.currentFile(1:end-7)],handles.trialsProjection);
%assignin('base', 'trials', handles.Trials);








% --- Executes on selection change in ProjectionItem.
function ProjectionItem_Callback(hObject, eventdata, handles)
% hObject    handle to ProjectionItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ProjectionItem contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ProjectionItem


% --- Executes during object creation, after setting all properties.
function ProjectionItem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ProjectionItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uipanel6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over displayMaxProjectionRadioButton.
function displayMaxProjectionRadioButton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to displayMaxProjectionRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Value',1);
handles.displayChannel='ProjectAcrossTrials';
guidata(handles.figure1,handles) ;


% --- Executes on button press in displayMaxProjectionRadioButton.
function displayMaxProjectionRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to displayMaxProjectionRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayMaxProjectionRadioButton
set(hObject,'Value',1);
handles.displayChannel='ProjectAcrossTrials';
guidata(handles.figure1,handles) ;


% --- Executes on selection change in ProjecAcrossTrialsSecondStep.
function ProjecAcrossTrialsSecondStep_Callback(hObject, eventdata, handles)
% hObject    handle to ProjecAcrossTrialsSecondStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ProjecAcrossTrialsSecondStep contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ProjecAcrossTrialsSecondStep


% --- Executes during object creation, after setting all properties.
function ProjecAcrossTrialsSecondStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ProjecAcrossTrialsSecondStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TrimRows_Callback(hObject, eventdata, handles)
% hObject    handle to TrimRows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TrimRows as text
%        str2double(get(hObject,'String')) returns contents of TrimRows as a double


% --- Executes during object creation, after setting all properties.
function TrimRows_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrimRows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TrimColumns_Callback(hObject, eventdata, handles)
% hObject    handle to TrimColumns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TrimColumns as text
%        str2double(get(hObject,'String')) returns contents of TrimColumns as a double


% --- Executes during object creation, after setting all properties.
function TrimColumns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrimColumns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderTrials_Callback(hObject, eventdata, handles)
% hObject    handle to sliderTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

try
    delete(handles.imageHandle)
catch
end


if ~isfield(handles,'trialsProjection')
    disp('Calculate Projection across trials first. Switching to Green Single Frame')
    activeChannel=handles.channels.G;
else
    activeChannel=handles.trialsProjection;
end

endingString=get(handles.endingTrialBox,'String');


trial=round(get(hObject,'Value'));
im=showImage(activeChannel(:,:,trial),handles);
handles.imageHandle=im;
uistack(handles.imageHandle,'bottom')

set(handles.TrialCounter,'String',[num2str(trial) '/' endingString ]);
guidata(handles.figure1,handles) ;

% --- Executes during object creation, after setting all properties.
function sliderTrials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject,'Min',1,'Max',5,'Value',1,'SliderStep',[ 1/5 10/5]);



function ComponentNum_Callback(hObject, eventdata, handles)
% hObject    handle to ComponentNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ComponentNum as text
%        str2double(get(hObject,'String')) returns contents of ComponentNum as a double


% --- Executes during object creation, after setting all properties.
function ComponentNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ComponentNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AutoROIsICA.
function AutoROIsICA_Callback(hObject, eventdata, handles)
% hObject    handle to AutoROIsICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try;
    delete(handles.roiPlotHandle(:),handles.roiTextHandle(:))
catch
    handles.roiPlotHandle=[];
    handles.roiTextHandle=[];
end
handles.roiPlotHandle=[];
handles.roiTextHandle=[];
handles.roipolyCoord=[];


ROIminSize=str2double(get(handles.ROIminSize,'String'));
contourWidth=str2double(get(handles.contourWidth,'String'));
trimRows=str2double(get(handles.TrimRows,'String'));
trimColumns=str2double(get(handles.TrimColumns,'String'));

SVDCompNumMax=str2double(get(handles.ComponentNum,'String'));
ICompNum=str2double(get(handles.IComponentNum,'String'));

if ~isfield(handles,'PCAData')
    
    disp('Loading PCA Data....');
    handles.PCAData=load( fullfile(handles.currentDirectory,['PCAData_' handles.currentFile(1:end-7)]),'Data');
    handles.U=load( fullfile(handles.currentDirectory,['PCAData_' handles.currentFile(1:end-7)]),'U');
    handles.S=load( fullfile(handles.currentDirectory,['PCAData_' handles.currentFile(1:end-7)]),'S');
    handles.V=load( fullfile(handles.currentDirectory,['PCAData_' handles.currentFile(1:end-7)]),'V');
end
[componentMasks L]=ICABasedROIAnalysis(handles.PCAData.Data,handles.U.U,handles.S.S,handles.V.V,handles.numberOfRows,handles.numberOfColumns,SVDCompNumMax,ICompNum);


trimLi=L;
trimLi(1:trimRows,:)=0;
trimLi(handles.numberOfRows-trimRows:end,:)=0;
trimLi(:,1:trimColumns)=0;
trimLi(:,handles.numberOfColumns-trimColumns:end)=0;

se90 = strel('line', contourWidth, 90);
se0 = strel('line', contourWidth, 0);
ROIindex=unique(trimLi(trimLi>0));
hold on
trimL=trimLi;
for i=1:length(ROIindex);
    if numel(find(trimL==ROIindex(i)))<ROIminSize
        trimL(trimL==ROIindex(i))=0;
    end
end

ROIindex=unique(trimL(trimL>0));
handles.roiNumber=length(ROIindex);
for i=1:length(ROIindex);
    I0=zeros(size(L));
    I0=trimL==ROIindex(i);
    I0dil=imdilate(I0,[ se90 se0] );
    ROIcontour=contourc(I0dil,1);
    handles.roipolyCoord{i}=(ROIcontour(:,2:end))';
    handles.roiPlotHandle(i)=plot(handles.roipolyCoord{i}(:,1),handles.roipolyCoord{i}(:,2),'Color',handles.colorlist(i,:),'LineWidth',2);
    handles.roiTextHandle(i)=text(mean(handles.roipolyCoord{i}(:,1))-2*length(num2str(i)), mean(handles.roipolyCoord{i}(:,2)), num2str(i),'Color',handles.colorlist(i,:),'FontWeight','bold');
end


set(handles.listboxROIs,'String',num2str([1:handles.roiNumber]'),...
    'Value',1);
uistack(handles.imageHandle,'bottom')
guidata(handles.figure1,handles) ;

uistack(handles.imageHandle,'bottom')
guidata(handles.figure1,handles) ;


function PCAMaxTrials_Callback(hObject, eventdata, handles)
% hObject    handle to PCAMaxTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PCAMaxTrials as text
%        str2double(get(hObject,'String')) returns contents of PCAMaxTrials as a double


% --- Executes during object creation, after setting all properties.
function PCAMaxTrials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PCAMaxTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IComponentNum_Callback(hObject, eventdata, handles)
% hObject    handle to IComponentNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IComponentNum as text
%        str2double(get(hObject,'String')) returns contents of IComponentNum as a double


% --- Executes during object creation, after setting all properties.
function IComponentNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IComponentNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PCAstep.
function PCAstep_Callback(hObject, eventdata, handles)
% hObject    handle to PCAstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PCAMaxTrials=str2double(get(handles.PCAMaxTrials,'String'));
[Data,U,S,V]=PrepareData(handles.currentFile(1:end-7),PCAMaxTrials,handles.numberOfRows,handles.numberOfColumns);

save( fullfile(handles.currentDirectory,['PCAData_' handles.currentFile(1:end-7)]),'Data','U','S','V','-v7.3');
if isfield(handles,'PCAData');
    rmfield(handles,[{'PCAData'};{'U'};{'S'};{'V'}])
end






% function ImClickCallback(hObject, eventdata,handles)
%
% if  get(handles.DonutAutoOn,'Value')  % only if donut check is on
%
%         lastROI=get(handles.listboxROIs,'String');
%         handles.roiNumber=str2num(lastROI(end,:))+1;
%         cellRadius=str2double(get(handles.donutCellRadius,'String'));
%         expandFactor=str2num(get(handles.donutExpandFactor,'String'));
% set(handles.listboxROIs,'String',num2str([1:handles.roiNumber]'),...
%     'Value',handles.roiNumber)
%     roiActive=handles.roiNumber;
%         pt = round(get(handles.axes1,'CurrentPoint'));
%         currentImage=get(handles.imageHandle,'CData');
%         I=round(pt(1,2));
%         J=round(pt(1,1));
% %       x=403;
% %       y=57;
%         [pixel_list,nuc_list,boundary_inner,boundary_outer]=donut_roi(currentImage,[I,J],cellRadius,expandFactor,0);
%
%         MaxString=get(handles.MaxValueText,'String');
%         ColormapMax=str2num(MaxString(2:end));
%         hax=handles.axes1;
%         hold(hax,'on');
%         m=zeros(size(currentImage));
%         m(pixel_list)=ColormapMax*2;
%         %him=imagesc(m,'AlphaData',m>1,'HitTest','off');
%         hi=plot(hax,boundary_inner(1,:),boundary_inner(2,:),'g');
%         ho=plot(hax,boundary_outer(1,:),boundary_outer(2,:),'g');
%         %set(handles.axes1,'XDir','Normal','YDir','Normal');
%         htext=text(mean(boundary_outer(1,:))-2*(length(num2str(roiActive))), mean(boundary_outer(2,:)), num2str(roiActive),'Color',handles.colorlist(roiActive,:),'FontWeight','bold','Parent', handles.axes1);
%
%         handles.roiTextHandle(roiActive)=htext;
%         handles.roiPlotHandle{handles.roiNumber}=[ho hi];
%         handles.ROIPixelList{handles.roiNumber}.outer=pixel_list;
%         handles.ROIPixelList{handles.roiNumber}.nuclear=nuc_list;
%         handles.roipolyCoord{roiActive}.outer=[ boundary_outer(1,:); boundary_outer(2,:)]';
%         handles.roipolyCoord{roiActive}.inner=[ boundary_inner(1,:); boundary_inner(2,:)]';
%         guidata(handles.figure1,handles) ;
%         %set(handles.imageHandle,'ButtonDownFcn',{@ImClickCallback, handles});
%      %  handles.ROIPixelList
% end




% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over listboxROIs.
function listboxROIs_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to listboxROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on listboxROIs and none of its controls.
function listboxROIs_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listboxROIs (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uitoggletool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.imageHandle,'ButtonDownFcn',{})



function donutCellRadius_Callback(hObject, eventdata, handles)
% hObject    handle to donutCellRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of donutCellRadius as text
%        str2double(get(hObject,'String')) returns contents of donutCellRadius as a double


% --- Executes during object creation, after setting all properties.
function donutCellRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to donutCellRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function donutExpandFactor_Callback(hObject, eventdata, handles)
% hObject    handle to donutExpandFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of donutExpandFactor as text
%        str2double(get(hObject,'String')) returns contents of donutExpandFactor as a double


% --- Executes during object creation, after setting all properties.
function donutExpandFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to donutExpandFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in setDonutROI.
function setDonutROI_Callback(hObject, eventdata, handles)
% hObject    handle to setDonutROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pt=ginput(1);

pt = round(get(handles.axes1,'CurrentPoint'));
% Changed by Tiago
%       lastROI=get(handles.listboxROIs,'String');
%       handles.roiNumber=str2num(lastROI(end,:))+1;
%         set(handles.listboxROIs,'String',num2str([1:handles.roiNumber]'),...
%         'Value',handles.roiNumber)
handles.roiNumber=get(handles.listboxROIs,'Value');

cellRadius=str2double(get(handles.donutCellRadius,'String'));
expandFactor=str2num(get(handles.donutExpandFactor,'String'));

roiActive=handles.roiNumber;


try
    delete(handles.roiPlotHandle{roiActive},handles.roiTextHandle(roiActive),handles.ROIPixelList{roiActive}.inner )
catch
end


% pt = round(get(handles.axes1,'CurrentPoint'));
currentImage=get(handles.imageHandle,'CData');
I=round(pt(1,2));
J=round(pt(1,1));
%       x=403;
%       y=57;
[pixel_list,nuc_list,boundary_inner,boundary_outer]=donut_roi(currentImage,[I,J],cellRadius,expandFactor,0);

MaxString=get(handles.MaxValueText,'String');
ColormapMax=str2num(MaxString(2:end));
hax=handles.axes1;
hold(hax,'on');
%m=zeros(size(currentImage));
%m(pixel_list)=ColormapMax*2;
%him=imagesc(m,'AlphaData',m>1,'HitTest','off');
hi=plot(hax,boundary_inner(1,:),boundary_inner(2,:),'g');
ho=plot(hax,boundary_outer(1,:),boundary_outer(2,:),'g');
%set(handles.axes1,'XDir','Normal','YDir','Normal');
htext=text(mean(boundary_outer(1,:)), mean(boundary_outer(2,:)), num2str(roiActive),'Color',handles.colorlist(roiActive,:),'FontWeight','bold','Parent', handles.axes1,'FontSize',6);

handles.roiTextHandle(roiActive)=htext;
handles.roiPlotHandle{handles.roiNumber}=[ho hi];
handles.ROIPixelList{handles.roiNumber}.outer=pixel_list;
handles.ROIPixelList{handles.roiNumber}.nuclear=nuc_list;
handles.roipolyCoord{roiActive}.outer=[ boundary_outer(1,:); boundary_outer(2,:)]';
handles.roipolyCoord{roiActive}.inner=[ boundary_inner(1,:); boundary_inner(2,:)]';
guidata(handles.figure1,handles) ;
%set(handles.imageHandle,'ButtonDownFcn',{@ImClickCallback, handles});
%  handles.ROIPixelList





function CCThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to CCThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CCThreshold as text
%        str2double(get(hObject,'String')) returns contents of CCThreshold as a double


% --- Executes during object creation, after setting all properties.
function CCThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CCThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CalculateROI_CrossCorr.
function CalculateROI_CrossCorr_Callback(hObject, eventdata, handles)
% hObject    handle to CalculateROI_CrossCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch handles.displayChannel
    case 'G'
        activeChannel=double(handles.channels.G);
    case 'R'
        activeChannel=double(handles.channels.R);
end
if ~isfield(handles,'ccimage');
    disp('Calculating Cross Correlations....Please Wait');
    [ccimage]=CrossCorrImage(activeChannel);
    
    handles.ccimage=ccimage;
end
%----

% try;
%     delete(handles.roiPlotHandle{:},handles.roiTextHandle(:))
% catch
%     handles.roiPlotHandle=[];
%     handles.roiTextHandle=[];
% end
% handles.roiPlotHandle=[];
% handles.roiTextHandle=[];
% handles.roipolyCoord=[];

threshold=str2double(get(handles.CCThreshold,'String'));
ROIminSize=str2double(get(handles.ROIminSize,'String'));
contourWidth=str2double(get(handles.contourWidth,'String'));
trimRows=str2double(get(handles.TrimRows,'String'));
trimColumns=str2double(get(handles.TrimColumns,'String'));

currentImage=handles.ccimage;
currentImage(1:trimRows,:)=0;
currentImage(handles.numberOfRows-trimRows:end,:)=0;
currentImage(:,1:trimColumns)=0;
currentImage(:,handles.numberOfColumns-trimColumns:end)=0;


BW=currentImage>threshold;
BWA=bwareaopen(BW,ROIminSize);
seD = strel('disk',contourWidth);
se90 = strel('line', contourWidth, 90);
se0 = strel('line', contourWidth, 0);


BWAF=(imdilate(imclearborder(imfill(BWA,'holes')),[ se90 se0]));
BWAFF=imclearborder(imfill(BWAF,'holes'),8);

BWfinal=imerode(BWAFF,seD);
BWfinal=bwareaopen(imerode(BWfinal,seD),ROIminSize);
%cc=bwconncomp(BWAF,8);
cc=bwlabel(BWfinal,8);

hold on
if isfield(handles,'roiTextHandle');
    currentNumberOfROIs=length(handles.roiTextHandle);
else
    
    currentNumberOfROIs=0;
end
for i=1:max(cc(:));
    roiIndex=i+currentNumberOfROIs;
    currentROIPix=double((cc==i));
    ROIcontour=contourc(currentROIPix,1);
    handles.roipolyCoord{roiIndex}.outer=(ROIcontour(:,2:end))';
    handles.roipolyCoord{roiIndex}.inner=[];
    handles.ROIPixelList{roiIndex}.outer=find(currentROIPix);
    handles.ROIPixelList{roiIndex}.inner=[];
    
    handles.roiPlotHandle{roiIndex}=plot(handles.roipolyCoord{roiIndex}.outer(:,1),handles.roipolyCoord{roiIndex}.outer(:,2),'Color',handles.colorlist(roiIndex,:),'LineWidth',2);
    handles.roiTextHandle(roiIndex)=text(mean(handles.roipolyCoord{roiIndex}.outer(:,1))-2*length(num2str(roiIndex)), mean(handles.roipolyCoord{roiIndex}.outer(:,2)), num2str(roiIndex),'Color',handles.colorlist(roiIndex,:),'FontWeight','bold');
end
handles.roiNumber=max(cc(:))+currentNumberOfROIs;

set(handles.listboxROIs,'String',num2str([1:handles.roiNumber]'),...
    'Value',1);
uistack(handles.imageHandle,'bottom')
guidata(handles.figure1,handles) ;

uistack(handles.imageHandle,'bottom')
guidata(handles.figure1,handles) ;


function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton36.
function pushbutton36_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AutoROIs.
function pushbutton37_Callback(hObject, eventdata, handles)
% hObject    handle to AutoROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in PlaneNumber.
function PlaneNumber_Callback(hObject, eventdata, handles)
% hObject    handle to PlaneNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlaneNumber contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlaneNumber
availablePlanes=get(handles.PlaneNumber,'String');
handles.planeOfInterest=str2num(availablePlanes(get(handles.PlaneNumber,'Value')));
guidata(handles.figure1,handles) ;


% --- Executes during object creation, after setting all properties.
function PlaneNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlaneNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in UseGPU.
function UseGPU_Callback(hObject, eventdata, handles)
% hObject    handle to UseGPU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UseGPU


% --- Executes on button press in LoadProjection.
function LoadProjection_Callback(hObject, eventdata, handles)
% hObject    handle to LoadProjection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'trialsProjection')
    handles.trialsProjection=[];
    
end

[f p]=uigetfile('Pick saved Projaction Across Trials');
l=load([p f]);
fnames=fieldnames(l);
handles.trialsProjection=l.(fnames{1});

string2=get(handles.ProjecAcrossTrialsSecondStep,'String');
projectionSecondStep=string2(get(handles.ProjecAcrossTrialsSecondStep,'Value'),:);


axes(handles.axes1);

delete(handles.imageHandle),
switch projectionSecondStep{1}
    case 'Max'
        im=showImage(max(handles.trialsProjection,[],3),handles);
    case 'Avg'
        im=showImage(nanmean(handles.trialsProjection,3),handles);
end
colormap(gray);
handles.imageHandle=im;
uistack(handles.imageHandle,'bottom')


%set sliderTrials
endingTrialNum=size(handles.trialsProjection,3);
trial=round(get(handles.sliderTrials,'Value'));
set(handles.sliderTrials,'Min',1,'Max',endingTrialNum,'Value',1,'SliderStep',[ 1/endingTrialNum 10/endingTrialNum]);
set(handles.TrialCounter,'String',[num2str(trial) '/' num2str(endingTrialNum) ]);
%handles.currentFile=filename;
guidata(handles.figure1,handles) ;



% --- Executes on button press in ROIActiveCheckbox.
function ROIActiveCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to ROIActiveCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ROIActiveCheckbox


% --- Executes on button press in CalculateNeuropil.
function CalculateNeuropil_Callback(hObject, eventdata, handles)
% hObject    handle to CalculateNeuropil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CalculateNeuropil


% --- Executes on button press in PlotNeuropilContour.
function PlotNeuropilContour_Callback(hObject, eventdata, handles)
% hObject    handle to PlotNeuropilContour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roiActive=get(handles.listboxROIs,'Value');

if isfield(handles,'currentROINeuropilContour')&&(handles.currentROINeuropilContour==1)
    delete(handles.currentROINeuropilPlotHandle);
    handles.currentROINeuropilContour=0;
    
else
    R=handles.ROIContourMaskAllROIs(:,:,roiActive);
    %R(handles.ROIPixelList{roiActive}.outer)=1;
    M=mask2poly(R,'MINDIST');
    goodPoints=M(:,1)>0;
    MGood=[M(goodPoints,1) M(goodPoints,2)];
    h=plot(MGood(1:end,1),MGood(1:end,2),'r');
    handles.currentROINeuropilContour=1;
    handles.currentROINeuropilPlotHandle=h;
end
guidata(handles.figure1,handles) ;

% --- Executes on button press in setCircleROI.
function setCircleROI_Callback(hObject, eventdata, handles)
% hObject    handle to setCircleROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%created by Marina 09.01.2015

set(handles.uitoggletool1,'State','off');
set(handles.uitoggletool3,'State','off');

pt=ginput(1);
handles.roiNumber=get(handles.listboxROIs,'Value');
roiActive=handles.roiNumber;
roiRadius=str2double(get(handles.donutCellRadius,'String'));
pt = round(get(handles.axes1,'CurrentPoint'));%current mouse point
cenx = pt(1,1); % center of circle x
ceny = pt(1,2); %center of circle y

%determine coordinates of circle to plot
x = roiRadius*cos(0:.1:2*pi)+cenx;
y = roiRadius*sin(0:.1:2*pi)+ceny;

%determine list of pixels inside circle
[indx, indy] = meshgrid(1:handles.numberOfRows,1:handles.numberOfColumns);
pixellist = find(((indx-cenx).^2+(indy-ceny).^2)<roiRadius^2); %assumes square imaging plane (equal # of rows and columns)


%plot newest ROI
hax=handles.axes1;
hold(hax,'on');
ho=plot(hax,x,y,'Color',handles.colorlist(roiActive,:));
htext=text(mean(x), mean(y), num2str(roiActive),'Color',handles.colorlist(roiActive,:),'FontWeight','bold','Parent', handles.axes1,'FontSize',6);
handles.roiTextHandle(roiActive)=htext;
handles.roiPlotHandle{handles.roiNumber}=[ho];
handles.ROIPixelList{handles.roiNumber}.outer=pixellist;
handles.ROIPixelList{handles.roiNumber}.nuclear=[];
handles.roipolyCoord{roiActive}.outer=[x; y]';
handles.roipolyCoord{roiActive}.inner=[];
guidata(handles.figure1,handles) ;


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

%created by Marina 24.05.2016
%same function as AddNewROI callback
if get(handles.figure1,'currentcharacter') == 'a'
    AddNewROI_Callback(hObject, eventdata, handles)
elseif get(handles.figure1,'currentcharacter') == 'c'
    setCircleROI_Callback(hObject, eventdata, handles)
elseif get(handles.figure1,'currentcharacter') == 'd' %added by gabi ag 2021 
    setDonutROI_Callback(hObject, eventdata, handles)
end

function [H,data] = load_SI2019a_EP(p,f)
[~, data]=opentif_SI5a(fullfile(p,f));
fprintf('ATTN: Couldn''t read header, extracting header from alternate file\n\n')
%opens second file to extract header
[f2,p2]=uigetfile('*.tif','Choose file with intact header');

% insert all the necessary header info here
[header,~,~,~] = scanimage.util.opentif(fullfile(p2,f2));
H.scanZoomFactor = header.SI.hRoiManager.scanZoomFactor;
H.scanLinePeriod = header.SI.hRoiManager.linePeriod;
H.scanFrameRate = header.SI.hRoiManager.scanFrameRate;
H.beamPowers = header.SI.hBeams.powers;
H.fastZEnable = header.SI.hFastZ.enable;
H.stackZStepSize = header.SI.hStackManager.stackZStepSize;
H.stackNumSlices = header.SI.hStackManager.numSlices;
H.channelsActive = header.SI.hChannels.channelsActive;
H.date = sprintf('%d%02d%02d',header.epoch{1}(1),header.epoch{1}(2),header.epoch{1}(3));
H.SI=header.SI;

function [H,data] = load_SI2019a_raw(p,f)
addpath(genpath('C:\Users\USER\Documents\MATLAB\Gabriela\FS_registration\SI2019bR1'));
[header,data,~] = scanimage.util.opentif(fullfile(p,f));

% insert all the necessary header info here
H.SI=header.SI;
H.scanZoomFactor = header.SI.hRoiManager.scanZoomFactor;
H.scanLinePeriod = header.SI.hRoiManager.linePeriod;
H.scanFrameRate = header.SI.hRoiManager.scanFrameRate;
H.beamPowers = header.SI.hBeams.powers;
H.fastZEnable = header.SI.hFastZ.enable;
H.stackZStepSize = header.SI.hStackManager.stackZStepSize;
if header.SI.hFastZ.enable %added Aug 2021 for neuron silencing data gabi
    H.stackNumSlices = header.SI.hStackManager.numSlices;
else
    H.stackNumSlices = 1;
end
H.channelsActive = header.SI.hChannels.channelSave;
H.date = sprintf('%d%02d%02d',header.epoch{1}(1),header.epoch{1}(2),header.epoch{1}(3));



% --- Executes on button press in bidicorrect_flag.
function bidicorrect_flag_Callback(hObject, eventdata, handles)
% hObject    handle to bidicorrect_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bidicorrect_flag
handles.bidi_correct_offset_flag = get(hObject,'Value');
guidata(handles.figure1,handles)


function bidiOffset_Callback(hObject, eventdata, handles)
% hObject    handle to bidiOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bidiOffset as text
%        str2double(get(hObject,'String')) returns contents of bidiOffset as a double
handles.bidiOffset = str2double(get(hObject,'String'));
guidata(handles.figure1,handles) 

% --- Executes during object creation, after setting all properties.
function bidiOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bidiOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function corr_image=biDiOffsetCorrection(off_image,shift_n)

corr_image=off_image;
if length(size(off_image))==2
    even_rows = corr_image(2:2:end,:);
    even_rows = circshift(even_rows,shift_n,2);
    corr_image(2:2:end,:)=even_rows;
elseif length(size(off_image))==3
    even_rows = corr_image(2:2:end,:,:);
    even_rows = circshift(even_rows,shift_n,2);
    corr_image(2:2:end,:,:)=even_rows;
else
    warning('Something is wrong...');
end


function im = showImage(x,handles)
L = get(gca,{'xlim','ylim'});  % Get axes limits to keep zoom

if handles.bidi_correct_offset_flag && handles.bidiOffset ~= 0
    corr_image=biDiOffsetCorrection(x,handles.bidiOffset);
    im = imagesc(corr_image,[handles.currentColorMapMinValue handles.currentColorMapMaxValue]);
else
    im = imagesc(x,[handles.currentColorMapMinValue handles.currentColorMapMaxValue]);
end
zoom(gcf,'reset')
%set(gca,{'xlim','ylim'},L) %reset to axes limits before plotting

