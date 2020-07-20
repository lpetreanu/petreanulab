nTrials = length(FAcrossTrials_reg3f_p1_c1_CMloop1_L55p__00.data);
nROIs = size(FAcrossTrials_reg3f_p1_c1_CMloop1_L55p__00.data{1},1);
nRep = length(FAcrossTrials_reg3f_p1_c1_CMloop1_L55p__00.data)*2;

for i = 1:nTrials
nFramesPerTrials(i) = min(size(FAcrossTrials_reg3f_p1_c1_CMloop1_L55p__00.data{i},2));
end
nFrames_min = min(nFramesPerTrials);
nFrames_max = max(nFramesPerTrials);

nFrames_min = nFrames_min - mod(nFrames_min,4);

data = cell(1,nTrials);
for i = 1:nTrials
    data{i} = FAcrossTrials_reg3f_p1_c1_CMloop1_L55p__00.data{i}(1:nROIs,1:nFrames_min);
end
data_concatenated = cat(2,data{1,:});

data_GrandMean = mean(data_concatenated,2);
data_sub = data_concatenated - data_GrandMean;

data_stacked = reshape(data_sub,nROIs,nFrames_min/4,nRep*2);

for i = 1:nRep*2
    if mod(i,2) == 0
        data_photostim(:,:,i/2) = data_stacked(:,:,i);
    else
        data_ctrl(:,:,floor(i/2)+1) = data_stacked(:,:,i);
    end
end

data_photostim_mean = nanmean(data_photostim,3);
data_photostim_sem = std(data_photostim,0,3)/sqrt(size(data_photostim,3));
Grand_mean = mean(data_photostim_mean,1);
for i = 1:nROIs
    figure; hold on
    plot(data_photostim_mean(i,:),'k')
    plot(data_photostim_mean(i,:)+data_photostim_sem(i,:),'Color',[0.5 0.5 0.5]);
    plot(data_photostim_mean(i,:)-data_photostim_sem(i,:),'Color',[0.5 0.5 0.5]);
    plot(Grand_mean,'Color','r')
end
% data_photostim_conc = reshape(data_photostim,nROIs,nFrames_min/4*nRep);