function PlotROIsMasks_2(data, rescale_factor,varargin)
global j
%% load data file to get the stat for the iscell cells


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rescale = 1; % if z-stack and functional recording have been acquired with different resolution
% rescale_factor = 1/2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if exist('path','var')
%     path = uigetdir(path,'directory to save ROIsMasks image');
if isempty(varargin)
    path = uigetdir('E:\TempData\CMloop33\20200603\pos3_zStack\Processed\','directory to save ROIsMasks image');
else
    path = varargin{1};
end
target = data.s2p{j}.meanImg;


ROIsMasks = zeros(size(target,1),size(target,2));
ROI_masks = figure; hold on
imagesc(target); colormap('gray')
for  i = 1:data.Info.nROIs(j)
    if data.s2p{j}.iscell(i)
        x = double(data.s2p{j}.stat{1, i}.xpix);
        y = double(data.s2p{j}.stat{1, i}.ypix);
                plot(x,y,'r');
        text(data.s2p{j}.stat{1, i}.med(2),data.s2p{j}.stat{1, i}.med(1),num2str(i),...
            'Color',[1 1 1],'HorizontalAlignment','center');
        for ii = 1:size(x,2)
            ROIsMasks(x(ii),y(ii)) = 1;
        end
%         for ii = 1:size(y,2)
%             test(y(ii)) = 1;
%         end
    end
end
axis equal
xlim([0 size(data.s2p{j}.meanImg,1)]); ylim([0 size(data.s2p{j}.meanImg,2)]);
axis ij

ROIsMasks = ROIsMasks';
figure; imagesc(ROIsMasks)

ROIsMasks = int16(ROIsMasks);
if rescale_factor ~= 1
ROIsMasks = imresize(ROIsMasks,rescale_factor);
end

oldfolder = cd(path);
tt=Tiff(['ROIsMasks_p' num2str(j) '.tif'],'w');
tagstruct.ImageLength = size(ROIsMasks,1);
tagstruct.ImageWidth = size(ROIsMasks,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.SampleFormat=Tiff.SampleFormat.Int;
tt.setTag(tagstruct);
tt.write(ROIsMasks);
tt.writeDirectory()
close(tt)
cd(oldfolder);

% % [composite_filename path] = uigetfile({'X:\camille.mazo\2P_processed\*.*'},'Load composite file');
% % composite = load([path composite_filename]);
% % composite(:,:,2) = test;
% 
% composite = source_reg_crop;
% composite(:,:,2) = source_c2_reg_crop;
% composite(:,:,3) = target;
% composite(:,:,4) = test;
% 
% composite = int16(composite);
% tt=Tiff('test.tif','w');
% tagstruct.ImageLength = size(composite,1);
% tagstruct.ImageWidth = size(composite,2);
% tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
% tagstruct.BitsPerSample = 16;
% tagstruct.SamplesPerPixel = 1;
% tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
% tagstruct.SampleFormat=Tiff.SampleFormat.Int;
% 
% 
% for ii = 1:size(composite,3)
%     for i = 1:size(composite,4)
%     tt.setTag(tagstruct);
%     tt.write(composite(:,:,ii,i));
%     %             cd(saveFolder_c1);
%     tt.writeDirectory()
%     %             cd ..
% end
% end
% close(tt)
end