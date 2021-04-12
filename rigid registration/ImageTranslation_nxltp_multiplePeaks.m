function [new_img shift] = ImageTranslation_nxltp_multiplePeaks(src_img, target_img ,padding, crop,varargin)
% Take the output from dftregistration.m and translate the image
% acoordingly. Creating a larger image accomadating the translation, no
% pixel cut off.

% shift, a 2-by-nFrame array of net_row_shift and net_col_shift
% padding, [top bottom left right] rows or columns to be padded to the src_img
%          if [], the 4 values would be taken from "shift".
%          If padding = [0 0 0 0], then the new_image is of the same size
%          of the original image, and the pixels shifted outside the image
%          are cut off.
% varargin: save_flag, save_path, save_name,ImageDescription

% The info for row and column index of the original image frames is saved
% in the header as 4-by-nframe array, with the 4 entries in each column to
% be [first_row_index; last_row_index; first_col_index; last_col_index], when save_flag is 1

% - NX 7/2009
%LP 4/2014 Modified to write using Tiff objects and support for Int16
%LP 8/2014
%M

if length(varargin)>=1
    save_flag = varargin{1};
else
    save_flag = 0;
end
if length(varargin)>=2
    save_path = varargin{2};
else
    save_path = pwd;
end
if length(varargin)>=3
    save_name = varargin{3};
else
    save_name = 'dft_reg_corrected';
end
% register every frame with dft reg algorithm, and get the shift value
shiftsOptionsX=NaN(5,size(src_img,3));
shiftsOptionsY=shiftsOptionsX;
if ~isempty(crop)
    target_img_cropped = target_img(crop(1):crop(2),crop(3):crop(4));
    src_img_cropped = src_img(crop(1):crop(2),crop(3):crop(4),:);
end
for i=1:size(src_img,3)
    [shiftsOptionsX(:,i),shiftsOptionsY(:,i)]  = dftregistrationMultiplePeaksLeo(fft2(double(target_img_cropped)),fft2(double(src_img_cropped(:,:,i))),1);
end

treshX=20;%<--------------------CONFIG
treshY=20;%<--------------------CONFIG

shift(1,:)=shiftsOptionsX(1,:);
shift(2,:)=shiftsOptionsY(1,:);
medianX=median(shift(1,:));
medianY=median(shift(2,:));
%edited by marina 18.11.2015 to add more options in case of big jumps
for i=2:5 %tries 5 options to avoid big jumps
    framesWrong=abs(shift(1,:)-medianX)>treshX|abs(shift(2,:)-medianY)>treshY;
    
    shift(1,framesWrong)=shiftsOptionsX(i,framesWrong);
    shift(2,framesWrong)=shiftsOptionsY(i,framesWrong);
end
% if there is still a big jump, take the shift value from the previous
% good frame (changed by Marina 09.10.2018)
framesWrong=find(abs(shift(1,:)-medianX)>treshX|abs(shift(2,:)-medianY)>treshY);
framesRight = setdiff(1:size(src_img,3),framesWrong);
for i = 1:length(framesWrong)
    s = find(framesRight<framesWrong(i),1,'last');
    if ~isempty(s)
        shift(1,framesWrong(i))=shift(1,framesRight(s));
        shift(2,framesWrong(i))=shift(2,framesRight(s));
    else
        shift(:,framesWrong(i)) = 0;
    end
end

% hack to not use suggested shift values and save non-registered images split into files by plane
%  shift = zeros(size(shift));

% median filter to remove outliers_______________________ LTP CONFIG VALUE !
%shift = round(medfilt1(squeeze(output(1,3:4,:))',2)');

% create a larger image according to the range of shift to be done
if isempty(padding)
    padding = get_im_padding(min(shift(1,:)),max(shift(1,:)),min(shift(2,:)),max(shift(2,:)));
end

class_str = class(src_img);
if length(size(src_img))==2
    new_img = zeros(size(src_img) + [padding(1)+padding(2), padding(3)+padding(4)], class_str)+mean(mean(src_img));
else
    new_img = zeros(size(src_img) + [padding(1)+padding(2), padding(3)+padding(4),0], class_str);
end
ind_row = padding(1)+1 : padding(1)+size(src_img,1);
ind_col = padding(3)+1 : padding(3)+size(src_img,2);


% get ImageDescription
if length(varargin) > 3
    im_describ = varargin{4};
else
    im_describ = '';
end
% If the shifts exceed the size if padded image, then cut off the exceeded
% pixels of the image.
%%

%tifStream
filename=[save_path filesep save_name];


%TiffObject
if length(size(src_img))>2

tt=Tiff(filename,'w');
tagstruct.ImageLength = size(src_img,1);
tagstruct.ImageWidth = size(src_img,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
% tagstruct.ImageDescription=im_describ.getTag('ImageDescription');
tagstruct.SampleFormat=Tiff.SampleFormat.Int;
%tagstruct.SubIFD=size(src_img,3);
tt.setTag(tagstruct);
end
%%


for i = 1:size(src_img,3)
    new_row = ind_row+shift(1,i);
    new_col = ind_col+shift(2,i);
    if new_row(1)<=0,
        top_cut = 1-new_row(1);
        new_row = 1:new_row(end);
    else
        top_cut = 0;
    end
    if new_row(end)>size(new_img,1)
        bottom_cut = new_row(end)-size(new_img,1);
        new_row = new_row(1):size(new_img,1);
    else
        bottom_cut =0;
    end
    if new_col(1)<=0
        left_cut = 1-new_col(1);
        new_col=1:new_col(end);
    else
        left_cut = 0;
    end
    if new_col(end)>size(new_img,2)
        right_cut = new_col(end)-size(new_img,2);
        new_col = new_col(1):size(new_img,2);
    else
        right_cut = 0;
    end
    
    new_img(new_row, new_col, i) = src_img((1+top_cut:end-bottom_cut),(1+left_cut:end-right_cut),i);
    
    
    if save_flag ==1
        if length(size(src_img))==2
            imwrite(uint16(new_img),filename)

        else
        header = im_describ;
        charFrameHeader = most.util.structOrObj2Assignments(header,'header');
        
        
        tagstruct.ImageLength = size(src_img,1);
        tagstruct.ImageWidth = size(src_img,2);
        % tagstruct.TileLength = size(src_img,1);
        % tagstruct.TileWidth = size(src_img,2);
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 16;
        tagstruct.SamplesPerPixel = 1;
        tagstruct.RowsPerStrip = 8;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.ImageDescription=charFrameHeader;
        tagstruct.SampleFormat=Tiff.SampleFormat.Int;
        tt.setTag(tagstruct);
        
        try
            tt.write(new_img(:,:,i));
            tt.writeDirectory()
        catch
            disp(['Unable to write frame # ' num2str(i)])
        end
        end
    end
end
if length(size(src_img))>2
    close(tt)
end

