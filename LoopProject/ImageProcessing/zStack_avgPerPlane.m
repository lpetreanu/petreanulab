function path = zStack_avgPerPlane(channel2, registration)
% load the raw z-stack file
% registration done on channel 1
% 
% channel2 = 1;
% registration = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[zStack_file, path] = uigetfile({'E:\TempData\*.*'},'Load a z-stack');
[~,stack]=opentif([path zStack_file]);
% 'E:\TempData\CMloop6\20200121_day8\CMloop6_da8_zstack\CMloop6_day8_L2_zStack_00001_00001.tif');
size(stack); %[x y c volumes frames optical_slices]
data_c2 = squeeze(stack(:,:,2,:,:,:)); 
data_c1 = squeeze(stack(:,:,1,:,:,:));

data_c2_avg = int16(squeeze(mean(data_c2,3)));
data_c1_avg = int16(squeeze(mean(data_c1,3)));

filename_c1 = [path zStack_file '_c1_avgPerPlane.tif'];
filename_c2 = [path zStack_file '_c2_avgPerPlane.tif'];

tt=Tiff(filename_c1,'w');
tagstruct.ImageLength = size(data_c1_avg,1);
tagstruct.ImageWidth = size(data_c1_avg,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.SampleFormat=Tiff.SampleFormat.Int;
if channel2
    tt2=Tiff(filename_c2,'w');
end
for ii = 1:size(data_c1_avg,3)
    tt.setTag(tagstruct);
    tt.write(data_c1_avg(:,:,ii));
    %             cd(saveFolder_c1);
    tt.writeDirectory()
    %             cd ..
end
close(tt)
if channel2
    for ii = 1:size(data_c2_avg,3)
        tt2.setTag(tagstruct);
        tt2.write(data_c2_avg(:,:,ii));
        %             cd(saveFolder_c2);
        tt2.writeDirectory()
        %             cd ..
    end
    close(tt2)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if registration
    for i = 1:size(data_c1_avg,3)
        fft_target = fft2(squeeze(data_c1_avg(:,:,i)));
        for ii = 1:size(data_c1,3)
            fft_source = fft2(squeeze(data_c1(:,:,ii,i)));
            [output ~] = dftregistration(fft_target,fft_source,1);
            x_shift(i,ii) = output(3);
            y_shift(i,ii) = output(4);
            error(i,ii) = output(1);
        end
    end
%     if channel2
%         for i = 1:size(data_c2_avg,3)
%             fft_target = fft2(squeeze(data_c2_avg(:,:,i)));
%             for ii = 1:size(data_c2,3)
%                 fft_source = fft2(squeeze(data_c2(:,:,ii,i)));
%                 [output ~] = dftregistration(fft_target,fft_source,1);
%                 x_shift(i,ii) = output(3);
%                 y_shift(i,ii) = output(4);
%                 error(i,ii) = output(1);
%             end
%         end
%     end


max_shift = max(max(max(abs(x_shift))),max(max(abs(y_shift))));
data_c1_reg = NaN(size(data_c1,1)+2*max_shift,size(data_c1,2)+2*max_shift,size(data_c1,3),size(data_c1,4));

for i = 1:size(data_c1,4)
    for ii = 1:size(data_c1,3)
%         if error(i,ii) < 0.3
            data_c1_reg(max_shift+x_shift(i,ii):size(data_c1,1)+max_shift+x_shift(i,ii)-1,max_shift+y_shift(i,ii):size(data_c1,2)+max_shift+y_shift(i,ii)-1,ii,i) = data_c1(:,:,ii,i);
%         else
%         end
    end
end
data_c1_reg_avg = squeeze(nanmean(data_c1_reg,3));

data_c1_reg_avg = int16(data_c1_reg_avg);
filename_c1_reg = [path zStack_file '_c1_avgPerPlane_reg.tif'];
tt=Tiff(filename_c1_reg,'w');
tagstruct.ImageLength = size(data_c1_reg_avg,1);
tagstruct.ImageWidth = size(data_c1_reg_avg,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.SampleFormat=Tiff.SampleFormat.Int;
for ii = 1:size(data_c1_reg_avg,3)
    tt.setTag(tagstruct);
    tt.write(data_c1_reg_avg(:,:,ii));
    %             cd(saveFolder_c1);
    tt.writeDirectory()
    %             cd ..
end
close(tt)
end