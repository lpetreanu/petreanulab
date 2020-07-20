% Simply give path to main position folder. Script is gonna import the
% different Fall values

% clearvars -except SessionData numberOfFrames saved_history IntrinsicCoord Ftraces_all;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_ROI = 1;

foldername = 'E:\TempData\CMloop33\20200605\pos4_stim'; % folder where all the files are located
% GenericFile_name = 'CMloop34_pos4*'; % 20190828_CMad47_pos1__00003
DataID = 'CMloop33_pos4';

path = [foldername '\suite2p\plane0'];

skipNofFrames = 1; % skip the extraction of number of frames from tiff files, coz it's an input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = datestr(now,'yyyymmdd');

%%
load([path '\Fall.mat'])

filetype = 'tif'; % type of files to be processed
saveFolder = [foldername '\analysis'];
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

oldfolder = cd(saveFolder);
diaryName = ['FtracesLog_' dt];
diary off; diary(diaryName)
cd(oldfolder)

saveName = ['Ftraces_' DataID '_' dt];
disp(saveFolder); disp(saveName);


% Y and X regstration offset
yOffset = ops.yoff1;
xOffset = ops.xoff1;
% list of files
fileList = ops.filelist;
        
if skipNofFrames==0
files = subdir(fullfile(foldername,[GenericFile_name,filetype]));   % list of filenames (will NOT search all subdirectories)
[H, Data] = opentif_SI2019b(files(1).name);
% H = opentif_SI2019b(files(1).name);

n2pFiles = length(files);
disp(['Found ',num2str(n2pFiles), ' tif files'])
numberOfFrames = zeros(n2pFiles,1);
f =  waitbar(0);
for i = 1:n2pFiles
    info = imfinfo(files(i).name);
    numberOfFrames(i) = numel(info)/(size(Data,3)*(size(Data,5)));
    waitbar(i/n2pFiles,f,[num2str(i/n2pFiles*100,2),'%'])
end
close(f)
else
    numberOfFrames = ops.frames_per_file;
end

disp(['number of frames per 2P files btwn ' num2str(max(numberOfFrames)) ' and ' num2str(min(numberOfFrames))])   

if all_ROI == 1
    Ftraces_all.data = F;
else
    Ftraces_all.data = F(is_cell(1:lastROI),:);
end

%% get header info from intact/splitted file
try
    Ftraces_all.header.frameFrequency=H.SI.hRoiManager.scanFrameRate;
catch
        indx = listdlg('PromptString','Sampling Frequency, per plane','SelectionMode','single',...
            'ListString',{'30','15','10','7.5','6'});
        Ftraces_all.header.frameFrequency =  30.4812/indx;
        disp(['Sampling Frequency = ' num2str(Ftraces_all.header.frameFrequency)]);
%    15.2406;
end

% Ftraces_all.header.SI5.stackNumSlices_original = H.SI.hStackManager.numSlices;
Ftraces_all.header.numberOfFrames = numberOfFrames;
% Ftraces_all.H = H;

Ftraces_all.Info.saveName = DataID;
Ftraces_all.Info.mainDir = foldername;
Ftraces_all.Info.fileList = fileList;

Ftraces_all.Fneu = Fneu;

Ftraces_all.s2p.iscell   = iscell(:,1);
Ftraces_all.s2p.meanImg  = ops.meanImg;
Ftraces_all.s2p.meanImgE = ops.meanImgE;
Ftraces_all.s2p.stat     = stat;
Ftraces_all.s2p.yOffset  = yOffset;
Ftraces_all.s2p.xOffset  = xOffset;

%% save 
fileName = [saveFolder,'\',saveName,'.mat'];
save(fileName,'Ftraces_all')
diary off