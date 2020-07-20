% [tiff1, folderPath1]=uigetfile('.tif','Select tif file 1');
% [tiff2, folderPath2]=uigetfile('.tif','Select tif file 2');
% captionText{1} = tiff1;
% captionText{2} = tiff2;

% t1 = Tiff([folderPath1 filesep tiff1],'r');
% t2 = Tiff([folderPath2 filesep tiff2],'r');

t1 = Tiff(['E:\TempData\CMloop19\z_stacks\pos4\CMloop19_pos4_zstack_ch830_00001_00001.tif_c1_avgPerPlane.tif_c1_reg.tif'],'r');
t2 = Tiff(['E:\TempData\CMloop19\z_stacks\pos4\CMloop19_pos4_zstack_ch830_00001_00001.tif_c2_avgPerPlane.tif_c2_reg.tif'],'r');


images{1} = read(t1);
images{2} = read(t2);


pixelInt_avg = NaN(length(data.s2p.stat),1);
pixelInt_max = NaN(length(data.s2p.stat),1);
bead_pos = false(length(data.s2p.stat),1);
bead_neg = false(length(data.s2p.stat),1);
for i = 1:length(data.s2p.stat)
    if data.s2p.iscell(i) == 1
        temp = images{1, 2}(data.s2p.stat{i}.xpix(1:2:end)/2,data.s2p.stat{i}.ypix(1:2:end)/2);
        temp = reshape(temp,size(temp,1)*size(temp,2),1);
        if any(data.beads_pos(:,1)==i) % i.e, bead+ cell
            bead_pos(i) = true;
            
        else
            bead_neg(i) = true;
        end
        pixelInt_avg(i,1) = sum(temp)/length(temp);
        pixelInt_max(i,1) = max(temp);
    end
end

figure;
subplot(2,1,1); hold on
histogram(pixelInt_avg(bead_neg))
histogram(pixelInt_avg(bead_pos))
title('Average pixel value within mask')
subplot(2,1,2); hold on
histogram(pixelInt_max(bead_neg))
histogram(pixelInt_max(bead_pos))
title('Highest pixel value within mask')
xlabel('Pixel Intensity')
%% to control for the masks
figure; imagesc(images{1, 2})
colormap('gray')
hold on
for i = 1: length(data.beads_pos  )
    plot(data.s2p.stat{data.beads_pos(i)}.xpix/2,data.s2p.stat{data.beads_pos(i)}.ypix/2,'r')
end
