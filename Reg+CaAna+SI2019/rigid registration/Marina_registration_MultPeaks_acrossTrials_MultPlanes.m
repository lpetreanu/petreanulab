function Marina_registration_MultPeaks_acrossTrials_MultPlanes(savefilename, channel, planesOfInterest, numberOfPlanes,ScanImage,byFrame)

%Registers a a series of movies against an specified image
% The list of files to register and teh target image are specified by uis.
%planesOfInterest is a vector ie: [1 2 3 4]
%number of PLanes, total number of Planes in experiment
%ScanImage - SI version
%byFrame: 1 - register frame by frame (1)  
%         0 - register trial  by trial
% Leopoldo Petreanu 2009
%%
%GF 2015
%SI5b, MP, mismatch

%GF 2016
%SI5a,MP,mismatch,MC
%%
%%
[f p] = uigetfile('*.tif','Select one of the files to be analyzed.');
cd(p)
[pathstr, name, ext] = fileparts([p f]);
d = dir([pathstr filesep '/*' ext]);

str = {d.name};
str = sortrows({d.name}');
[s,v] = listdlg('PromptString','Select files to analyze:', 'OKString', 'OK',...
    'SelectionMode','multiple',...
    'ListString', str, 'Name', 'Select a File');
names = str(s);
numTraces = size(names, 1);

for i=1:length(planesOfInterest)
    [targetfilename p]=uigetfile('*.tif*', ['Select Image Target File of PLANE# ' num2str(planesOfInterest(i))]);
    targets(:,:,i)=imread(targetfilename);
end

output=NaN(2,2000,numTraces);
[r, c, z]=size(targets);

% %frame-by-frame
MaxProjection=NaN(r,c,numTraces,z);

for i=1:numTraces
    fprintf('processing trial %d \n', i)
    fullsourcefilename = [pathstr filesep names{i}];
    
    [movAvgim,shift{i}]= Marina_registration_MultPeak_singleTrial_padding_multPlanes(targets, fullsourcefilename,savefilename,channel,planesOfInterest, numberOfPlanes,ScanImage,byFrame);
    MaxProjection(:,:,i,:)=movAvgim;
end

%frame-by-frame
%Declaring and formatting
%     MaxProjection=uint16(MaxProjection);

C=mat2cell(ones(numTraces,1)*2,repmat(1,1,numTraces),1);
trialsFrames=cellfun(@size,shift,C');
sizeToPad=max(trialsFrames);
ShiftPadded=NaN(2,sizeToPad,numTraces);


for i=1:length(planesOfInterest)
    %frame-by-frame
    %save MaxProjection of individual planes
    MaxProjectionPlane=squeeze(MaxProjection(:,:,:,i));
    save([pathstr filesep savefilename name(1:end-3) 'MaxProjection_Plane_' num2str(planesOfInterest(i))], 'MaxProjectionPlane');
    
    
    %Save XYSHift of individual planes
    if byFrame
        for j=1:numTraces;
            trialLenght=size(shift{j},2);
            ShiftPadded(:,1:trialLenght,j)=squeeze(shift{j}(:,:,i)); %[X Y],time,Trial
            save([pathstr filesep savefilename name(1:end-3) 'ShiftXY_Plane_' num2str(planesOfInterest(i))], 'ShiftPadded');
        end
    else
        
    end
end


figure(78);subplot(1,2,1);imagesc(squeeze(ShiftPadded(1,:,:))');colorbar;subplot(1,2,2);imagesc(squeeze(ShiftPadded(2,:,:))');colorbar;