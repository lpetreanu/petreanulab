function [new_img, shift] = ImageTranslation_nxltp_multiplePeaks_TrialbyTrial(src_img, target_img ,padding, SIver, crop, varargin)
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

trialProj = mean(src_img,3);
if ~isempty(crop)
    target_img = target_img(crop(1):crop(2),crop(3):crop(4));
    trialProj = trialProj(crop(1):crop(2),crop(3):crop(4));
end
% calculate amount of shift from the projection of each single trial
[shiftsOptionsX,shiftsOptionsY]  = dftregistrationMultiplePeaksLeo(fft2(double(target_img)),fft2(double(trialProj)),1);

treshX=20;%<--------------------CONFIG
treshY=20;%<--------------------CONFIG


%edited by marina 09.10.2018 to add more options in case of big jumps
i = 1;
while (abs(shiftsOptionsX(i))>treshX||abs(shiftsOptionsY(i))>treshY) && i<5 %tries 5 options to avoid big jumps
    i = i+1;
end
shift(1) = shiftsOptionsX(i);
shift(2) = shiftsOptionsY(i);

% hack to not use suggested shift values and save non-registered images split into files by plane
%  shift = zeros(size(shift));

% median filter to remove outliers_______________________ LTP CONFIG VALUE !
%shift = round(medfilt1(squeeze(output(1,3:4,:))',2)');

% create a larger image according to the range of shift to be done
if isempty(padding)
    padding = get_im_padding(shift(1),shift(1),shift(2),shift(2));
end

class_str = class(src_img);
if length(size(src_img))==2
    new_img = zeros(size(src_img) + [padding(1)+padding(2), padding(3)+padding(4)], class_str);
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

    %%


for i = 1:size(src_img,3)
    new_row = ind_row+shift(1);
    new_col = ind_col+shift(2);
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
%     if ~isequal(padding, [0 0 0 0])
%         % if the new image size is different than the original image then
%         % Add the index info of the original image to the header of the scanimgage file.
%         header_append_line1 = ['state.OriginalImageIndex.row1=' num2str(new_row(1))];
%         header_append_line2 = ['state.OriginalImageIndex.row_end=' num2str(new_row(end))];
%         header_append_line3 = ['state.OriginalImageIndex.col1=' num2str(new_col(1))];
%         header_append_line4 = ['state.OriginalImageIndex.col_end=' num2str(new_col(end))];
%         im_describ = [im_describ header_append_line1 char(13) header_append_line2 char(13) ...
%             header_append_line3 char(13) header_append_line4];
%     end
    

    if save_flag ==1
        % if i == 1,
        %     imwrite(new_img(:,:,i),[save_path filesep save_name],'tif','Compression','none','Description',im_describ,'WriteMode','overwrite');
        % else
        %imwrite(uint16(new_img(:,:,i)),[save_path filesep save_name],'tif','Compression','none','Description',im_describ,'WriteMode','append');
        % end;
        
        %% modify ScanImage Header
        switch SIver
            case 'SI5'
                if isfield(im_describ,'SI5'); %Versions before 2015a ( before August 2015)
                    % Fake DATA to fix defective openTif LTP
                    im_describ.nextFileMarkerTimestamps=im_describ.acqTriggerTimestamps;
                    im_describ.dcOverVoltage=im_describ.acqTriggerTimestamps;
                    
                    s=im_describ;
                    s.frameNumbers=im_describ.frameNumbers(i*plane);
                    s.frameTimestamps=im_describ.frameTimestamps(i*plane);
                    s.acqTriggerTimestamps=im_describ.acqTriggerTimestamps(i*plane);
                    s.nextFileMarkerTimestamps=im_describ.nextFileMarkerTimestamps(i*plane);
                    s.dcOverVoltage=im_describ.dcOverVoltage(i*plane);
                    
                    [snew, perm]=orderfields(s,{'frameNumbers','frameTimestamps','acqTriggerTimestamps','nextFileMarkerTimestamps','dcOverVoltage','acquisitionNumbers','frameNumberAcquisition','SI5'});
                    
                    s=snew;
                    
                    s.SI5.channelsSave=1;  %%<------------single channel for now. Assigned as GREEN)
                    s.SI5.stackNumSlices=1;
                    s.SI5.fastZEnable=1; % remove plane information
                    s.SI5.fastZNumVolumes=1;
                    
                    s=rmfield(s,{'acquisitionNumbers', 'frameNumberAcquisition'});
                    charFrameHeader=most.util.structOrObj2Assignments(s,'scanimage');
                    charFrameHeader=regexprep(charFrameHeader,'scanimage.frameNumbers','Frame Number');
                    charFrameHeader=regexprep(charFrameHeader,'scanimage.frameTimestamps','Frame Timestamp(s)');
                    charFrameHeader=regexprep(charFrameHeader,'scanimage.acqTriggerTimestamps','Acq Trigger Timestamp(s)');
                    charFrameHeader=regexprep(charFrameHeader,'scanimage.nextFileMarkerTimestamps','Next File Marker Timestamp(s)');
                    charFrameHeader=regexprep(charFrameHeader,'scanimage.dcOverVoltage','DC Overvoltage');
                else    % SI2015a August 2015, Installed Sept 2015
                    s=im_describ;
                    plane = 1;
                    s.frameNumbers=im_describ.frameNumbers(i*plane);
                    s.frameTimestamps=im_describ.frameTimestamps(i*plane);
                    s.acqTriggerTimestamps=im_describ.acqTriggerTimestamps(i*plane);
                    s.nextFileMarkerTimestamps=im_describ.nextFileMarkerTimestamps(i*plane);
                    s.endOfAcquisition=im_describ.endOfAcquisition(i*plane);
                    s.endOfAcquisitionMode=im_describ.endOfAcquisitionMode(i*plane);
                    s.dcOverVoltage=im_describ.dcOverVoltage(i*plane);
                    s.SI.hChannels.channelSave=1;  %%<------------single channel for now. Assigned as GREEN)
                    s.SI.stackNumSlices=1;
                    s.SI.fastZEnable=1; % remove plane information
                    s.SI.fastZNumVolumes=1;
                    
                    s=rmfield(s,{'acquisitionNumbers', 'frameNumberAcquisition'});
                    charFrameHeader=most.util.structOrObj2Assignments(s,'scanimage');
                    charFrameHeader=regexprep(charFrameHeader,'scanimage.frameNumbers','Frame Number');
                    charFrameHeader=regexprep(charFrameHeader,'scanimage.frameTimestamps','Frame Timestamp(s)');
                    charFrameHeader=regexprep(charFrameHeader,'scanimage.acqTriggerTimestamps','Acq Trigger Timestamp(s)');
                    charFrameHeader=regexprep(charFrameHeader,'scanimage.nextFileMarkerTimestamps','Next File Marker Timestamp(s)');
                    charFrameHeader=regexprep(charFrameHeader,'scanimage.dcOverVoltage','DC Overvoltage');
                end
            case 'SI4'
                header = im_describ;
                charFrameHeader = most.util.structOrObj2Assignments(header,'header');
        end
        %%
        
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
        
        
        
        tt.write(new_img(:,:,i));
        tt.writeDirectory()
        
    end
end
close(tt)


