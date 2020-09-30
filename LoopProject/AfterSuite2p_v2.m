%%%%%%%%%%% INSTEAD of saving data, Fneu, iscell, meanImg etc in different
%%%%%%%%%%% cells, save the whole thing in 2 different cells


% Simply give path to main position folder. Script is gonna import the
% different Fall values

% clear all
% clearvars -except SessionData numberOfFrames saved_history IntrinsicCoord Ftraces_all;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_ROI = 1; % I usually kee this at 1, selection is done in the next function

foldername = 'E:\TempData\CMloop44\20200907\pos1_stim'; % folder where all the files are located
% GenericFile_name = 'CMloop34_pos4*'; % 20190828_CMad47_pos1__00003
DataID = 'CMloop44_pos1';
skipNofFrames = 1; % to skip the extraction of number of frames from tiff files (using SI tiff reader), use s2p output instead

isFlyBack = 1; % Only matters if multiple planes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parametrize
s2p_folder = dir([foldername filesep 'suite2p']);
match = {}; counter = 1;
for i = 1:length(s2p_folder)
    if ~isempty(regexp(s2p_folder(i).name,'^plane'))
        match{end+1} = s2p_folder(i).name;
        index(counter) = i;
        counter = counter +1;
    end
end
if length(index)>1
    index = index(1:end-isFlyBack);
    plane = 'multiples';
else
    plane = 'single';
end

dt = datestr(now,'yyyymmdd');

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

disp(['number of planes: ', num2str(length(index))])
disp(saveFolder); disp(saveName);
disp(['All ROIs?: ', num2str(all_ROI)])

%% Chopp all the Fall file of each plane
for j = 1:length(index)
    path = [s2p_folder(index(j)).folder filesep s2p_folder(index(j)).name];
    disp(s2p_folder(index(j)).name)
    load([path '\Fall.mat'])
    
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
        Ftraces_all{j}.data = F;
    else
        Ftraces_all{j}.data = F(logical(iscell(:,1)),:);
    end
    
    Ftraces_all{j}.Fneu = Fneu;
    
    Ftraces_all{j}.s2p.iscell   = iscell(:,1);
    Ftraces_all{j}.s2p.meanImg  = ops.meanImg;
    Ftraces_all{j}.s2p.meanImgE = ops.meanImgE;
    Ftraces_all{j}.s2p.stat     = stat;
    Ftraces_all{j}.s2p.yOffset  = yOffset;
    Ftraces_all{j}.s2p.xOffset  = xOffset;
    
    
    try
        Ftraces_all{j}.header.frameFrequency=H.SI.hRoiManager.scanFrameRate;
    catch
        if j == 1
        indx = listdlg('PromptString','Sampling Frequency, per plane','SelectionMode','single',...
            'ListString',{'30','15','10','7.5','6'});
        frameFrequency =  30.4812/indx;
        disp(['Sampling Frequency = ' num2str(frameFrequency)]);
        end
    Ftraces_all{j}.header.frameFrequency = frameFrequency;
        
        %    15.2406;
    end
    % Ftraces_all.header.SI5.stackNumSlices_original = H.SI.hStackManager.numSlices;
    Ftraces_all{j}.header.numberOfFrames = numberOfFrames;
    % Ftraces_all.H = H;
    
    Ftraces_all{j}.Info.saveName = DataID;
    Ftraces_all{j}.Info.mainDir = foldername;
    Ftraces_all{j}.Info.fileList = fileList;
    Ftraces_all{j}.Info.plane = plane;
end


%% save
fileName = [saveFolder,'\',saveName,'.mat'];
save(fileName,'Ftraces_all')
diary off