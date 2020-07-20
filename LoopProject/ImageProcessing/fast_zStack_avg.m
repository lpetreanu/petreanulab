answer = inputdlg('How many Fast z-stack?');
k = str2num(answer{1,1});

dir = uigetdir('E:\TempData\','Select folder');
for i = 1:k
    [stack_fileName, path] = uigetfile([dir '\*.*'],'Load a fast z-stack');
    paths{i} = path;
    stack_fileNames{i} = stack_fileName;
end

data_c1_avg = NaN(size(data_c1_avg));
data_c2_avg = NaN(size(data_c2_avg));
for i = 1:k
    disp(stack_fileNames{i})
    [~,zStack]=opentif([paths{i} stack_fileNames{i}]);
    data_c1 = squeeze(zStack(:,:,1,:,:,:));
    data_c2 = squeeze(zStack(:,:,2,:,:,:));
    
    data_c1_avg(:,:,:,i) = mean(data_c1,4);
    data_c2_avg(:,:,:,i) = mean(data_c2,4);
end


%% save tiffs, separating plane 1 and 2
saveDir = [path 'Processed' filesep ];
if ~ exist(saveDir)
    mkdir(saveDir)
end

for i = 1:size(data_c1,3)
    dataToSave_c1 = int16(squeeze(data_c1_avg(:,:,i,:)));
    filename_c1 = [saveDir stack_fileName(1:end-4) '_c1_p' num2str(i) '.tif'];
    tt=Tiff(filename_c1,'w');
    tagstruct.ImageLength = size(dataToSave_c1,1);
    tagstruct.ImageWidth = size(dataToSave_c1,2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.SampleFormat=Tiff.SampleFormat.Int;
    for ii = 1:size(dataToSave_c1,3)
        tt.setTag(tagstruct);
        tt.write(dataToSave_c1(:,:,ii));
        tt.writeDirectory()
    end
    close(tt)
    
    dataToSave_c2 = int16(squeeze(data_c2_avg(:,:,i,:)));
    filename_c2 = [saveDir stack_fileName(1:end-4) '_c2_p' num2str(i) '.tif'];
    tt2=Tiff(filename_c2,'w');
    for ii = 1:size(dataToSave_c2,3)
        tt2.setTag(tagstruct);
        tt2.write(dataToSave_c2(:,:,ii));
        tt2.writeDirectory()
    end
    close(tt2)
end