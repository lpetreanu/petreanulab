function [maxZ,shift]=Marina_registration_MultPeak_singleTrial_padding_multPlanes(target_img, fullsourcefilename,savefilename,channel,planesOfInterest,numberOfPlanes,scanimage,byFrame)

% Image subpixel registration implamentation for Scan Image File
% Registers one channel only for now. 
% Based on dftregistration.m 
%Registered files are saved on the sane directory with prefix "reg_chan_1..."
% % Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).
% Leopoldo Petreanu 2009

%byFrame: 1 - register frame by frame (1)  
%         0 - register trial  by trial


% define the subpixel resoluton ( e.g. =10, subpixel rsolution =1/10)
subpixelFraction=1; %<---- Config
          
crop = [1 512 100 400]; %added 020818 to cut out lines created by monitor shutdown during scan
% crop = [1 200 50 200];

[pathstr, name, ext] = fileparts(fullsourcefilename);
cd(pathstr)
bail = 0; %option to bail on a particular trial if the file is corrupted
switch scanimage %changed for SI4 Gabi 270616
    case 'SI4'
        try
            [header,Stack] = scim_openTif(fullsourcefilename, 'channel',channel);
        catch
            disp(['Could not read file ' fullsourcefilename])
            bail = 1;
        end
    case 'SI5b'
        [~, header ,Stack] = Leo_opentif_sigleSliceSI5(fullsourcefilename,'channel',channel);
    case 'SI5a'
        [ header ,Stack] = opentif_SI5a(fullsourcefilename,'channel',channel,'slice',planesOfInterest);
end

if ~bail
    for i = 1:length(planesOfInterest)
        switch scanimage
            case 'SI4'
                G_Stack=squeeze(Stack(:,:,channel,planesOfInterest(i):numberOfPlanes:end));
                [r, c, z]=size(G_Stack);
                header.SI4.stackNumSlices=1;  % remove multPlaneData from Header since we're saving one plane per file
                header.SI4.fastZEnable=0;
                header.SI4.acqNumFrames=z;
                header.frameTags=header.frameTags(planesOfInterest(i):numberOfPlanes:end);
                header.SI4.channelsSave = 1;
            otherwise
                if channel==1
                    if size(Stack,3)==1;
                        G_Stack=squeeze(Stack(:,:,1,1,planesOfInterest(i),:));
                    else
                        G_Stack=Stack(:,:,planesOfInterest(i):numberOfPlanes:end);
                    end
                end
                if channel==2
                    if size(Stack,3)==1;
                        G_Stack=squeeze(Stack(:,:,2,1,planesOfInterest(i),:));
                    else
                        G_Stack=Stack(:,:,planesOfInterest(i):numberOfPlanes:end);
                    end
                end
        end
        
        savefilenamePlane=[num2str(i) filesep savefilename 'Plane_' num2str(i)];
        if byFrame
            % frame by frame
            [regImg, shiftTemp] = ImageTranslation_nxltp_multiplePeaks...
                (G_Stack, target_img(:,:,find(planesOfInterest==i)),[0 0 0 0],crop,1,pathstr, [savefilenamePlane  '_Chan' num2str(channel) '_' name '.tif'],header,planesOfInterest(i),1);
        else
            % trial by trial
            [regImg, shiftTemp] = ImageTranslation_nxltp_multiplePeaks_TrialbyTrial...
                (G_Stack, target_img(:,:,find(planesOfInterest==i)),[0 0 0 0],scanimage,crop,1,pathstr, [savefilenamePlane '_Chan' num2str(channel) '_' name '.tif'], header);
        end
        
        if i==1
            shift = shiftTemp;
        end
        if size(shiftTemp,2)==size(shift,2)
            shift(:,:,find(planesOfInterest==i)) = shiftTemp;
        else
            shift(:,:,find(planesOfInterest==i)) = nan;
        end
        movAvgim=im_mov_avg(regImg,5);
        maxZ(:,:,:,find(planesOfInterest==i))= max(movAvgim,[],3);
    end
else
    [r,c,z] = size(target_img);
    maxZ = nan(r,c,1,length(planesOfInterest));
    shift = nan(2,1,5);
end


