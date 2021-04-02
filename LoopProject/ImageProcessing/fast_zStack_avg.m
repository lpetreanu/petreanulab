answer = inputdlg('How many Fast z-stack?');
nVolumes = str2num(answer{1,1});

channel2 = 1;

dir = uigetdir('E:\TempData\','Select folder');
for j = 1:nVolumes
    [stack_fileName, path] = uigetfile([dir '\*.*'],'Load a fast z-stack');
    paths{j} = path;
    stack_fileNames{j} = stack_fileName;
end

% data_c1_avg = NaN(size(data_c1_avg));
% data_c2_avg = NaN(size(data_c2_avg));
for j = 1:nVolumes
    disp(stack_fileNames{j})
    [~,zStack]=opentif([paths{j} stack_fileNames{j}]);
    data_c1 = squeeze(zStack(:,:,1,:,:,:));
    data_c2 = squeeze(zStack(:,:,2,:,:,:));
    
    %% Self-register each plane
    % evaluate movement for plane1 only: assumes same motion between the 2 planes
    data_c1_p1Avg = squeeze(mean(data_c1(:,:,1,:),4));
    data_c1_p1 = squeeze(data_c1(:,:,1,:));
    nFrames = size(data_c1_p1,3);
    nPlanes = size(data_c1,3);
    
    for frame = 1:nFrames
        c1_fft = fft(data_c1_p1(:,:,frame));
        %    c1_mean(i) = mean2(c1_fft);
        %    c1_max(i) = max(max(c1_fft));
        %    plot(abs(mean(c1_fft,2)));
        crispy(:,frame,1) = abs(mean(c1_fft,1)); % max frequencies in x
        crispy(:,frame,2) = abs(mean(c1_fft,2)); % max frequencies in y
    end
    crispy_avg = mean(crispy,3); % average x and y fft powers
    
    for  frame =1:nFrames
        AUC(frame) = nansum(crispy_avg(350:end,frame));
    end
    figure; hold on
    histogram(AUC)
    threshold = 3000;
    plot([threshold threshold],ylim,'r')
    xlabel('Power of high frequencies')
    title('Cut off for planes selection')
    
    % Max projection excluding bad planes
    k = false(size(AUC,2),1);
    k(AUC>threshold) = true;
    figure;plot(k)
    xlabel('Plane #')
    title('Planes selected for Max Proj')
    % maxProj= max(data_c1_p1Avg(i,ii,k));
    maxProj= data_c1_p1Avg;
    
    % Register single frames to that Max Projection
    fft_target = fft2(maxProj);
    outputs = zeros(nFrames,4);
    for frame = 1:nFrames
        source = data_c1_p1(:,:,frame);
        fft_source = fft2(source);
        [output, ~] = dftregistration(fft_target,fft_source,1);
        outputs(frame,:) = output;
    end
    
    max_shift = max(max(max(abs(outputs(:,3)))), max(max(abs(outputs(:,4)))));
    source_reg_crop = NaN(size(source,1),size(source,2),nFrames,nPlanes);
    source_reg_crop_c2 = NaN(size(source,1),size(source,2),nFrames,nPlanes);
    
    for plane = 1:nPlanes
        for frame = 1:nFrames
            source_reg = NaN(size(source,1)+2*max_shift,size(source,2)+2*max_shift);
            source_reg(max_shift+1+outputs(frame,3):size(source,1)+max_shift+1+outputs(frame,3)-1,max_shift+1+outputs(frame,4):size(source,1)+max_shift+1+outputs(frame,4)-1) = data_c1(:,:,plane,frame);
            source_reg_crop(:,:,frame,plane) = source_reg(max_shift:max_shift+size(source,1)-1,max_shift:max_shift+size(source,2)-1);
            if channel2
                source_reg_c2 = NaN(size(source,1)+2*max_shift,size(source,2)+2*max_shift);
                source_reg_c2(max_shift+1+outputs(frame,3):size(source,1)+max_shift+1+outputs(frame,3)-1,max_shift+1+outputs(frame,4):size(source,1)+max_shift+1+outputs(frame,4)-1) = data_c2(:,:,plane,frame);
                source_reg_crop_c2(:,:,frame,plane) = source_reg_c2(max_shift:max_shift+size(source,1)-1,max_shift:max_shift+size(source,2)-1);
            end
        end
    end
    data_c1_avg(:,:,:,j) = squeeze(nanmean(source_reg_crop,3));
    data_c2_avg(:,:,:,j) = squeeze(nanmean(source_reg_crop_c2,3)); 
end

% %% Register the different volumes to themselves
% grandMean_p1 = mean(data_c1_avg(:,:,1,:),4);
% grandMean_p1 = data_c1_avg(:,:,1,3);
% 
% fft_target = fft2(grandMean_p1);
% for j = 1:nVolumes
%     fft_source = fft2(squeeze(data_c1_avg(:,:,1,j)));
%     [output, ~] = dftregistration(fft_target,fft_source,1);
%     outputs(j,:) = output;
% end

%% save tiffs, separating plane 1 and 2
saveDir = [path 'Processed' filesep ];
if ~ exist(saveDir)
    mkdir(saveDir)
end

for j = 1:nPlanes
    dataToSave_c1 = int16(squeeze(data_c1_avg(:,:,j,:)));
    filename_c1 = [saveDir stack_fileName(1:end-4) '_c1_p' num2str(j) '.tif'];
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
    
    dataToSave_c2 = int16(squeeze(data_c2_avg(:,:,j,:)));
    filename_c2 = [saveDir stack_fileName(1:end-4) '_c2_p' num2str(j) '.tif'];
    tt2=Tiff(filename_c2,'w');
    for ii = 1:size(dataToSave_c2,3)
        tt2.setTag(tagstruct);
        tt2.write(dataToSave_c2(:,:,ii));
        tt2.writeDirectory()
    end
    close(tt2)
end