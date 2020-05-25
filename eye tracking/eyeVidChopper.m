% for cutting eye video files down to size 
%Marina 17.11.14

%takes eye files saved by flycapture and which were triggered by scanimage
%frames. 3 eye files for every scanimage frame

%ask for eye file
[eyeMoviefn eyeMoviePath]=uigetfile('.avi','Select eye tracking movie file');
vidObj = VideoReader([eyeMoviePath filesep eyeMoviefn]);

%ask for trial that eye file is aligned with
firstTrial = str2num(input('What trial is this file aligned with? ','s'));

%ask for Facrosstrials file which has the number of frames of each file
[ffn fpath]=uigetfile('.mat','Select FAcrossTrials file');

%ask for base of eye file name
newfnbase = input('what shall I call the new files? ','s');

%gets number of trials and frames from FAcrossTrials file
f = load([fpath filesep ffn]);
fields = fieldnames(f);
f = f.(fields{1});
nt = size(f.data,2);
nfr = zeros(nt,1);
for i = firstTrial:nt
    nfr(i) = size(f.data{i},2);
end

%converts frame numbers from scanimage into frame numbers for eye file
indEnd = cumsum(nfr)*3;%last index of eyeframe for each trial
indStart = [1;indEnd+1];
indStart = indStart(1:end-1); %first index of eyeframe for each trial
lastTrial = find(indEnd<=vidObj.NumberOfFrames);
lastTrial = max(lastTrial);

%reads frames and saves to new file

for i = firstTrial:lastTrial
    
    if indStart(i) < indEnd(i) %only if there are frames in this trial
        %read desired frames
        fr = read(vidObj,[indStart(i) indEnd(i)]);
        fr = fr(:,:,1,:); %take only red channel
        fr = squeeze(fr);
        
        %analyze pupil
%         if i == firstTrial
%             figure,image(fr(:,:,1))
%             roiInd = round(getrect(gcf)); %returns indices of roi to analyze, [xmin ymin width height]
%         end
%         
%         roi = fr(roiInd(2):roiInd(2)+roiInd(4),roiInd(1):roiInd(1)+roiInd(3),:);
%         [pd,px,py] = analyzeThatPupil(roi);
        
        %write to new file
        newfn = [eyeMoviePath filesep newfnbase sprintf('%03d',i) '.avi'];
        newvid = VideoWriter(newfn);
        open(newvid);
        for j = 1:size(fr,3)
            writeVideo(newvid,squeeze(fr(:,:,j)))
        end
        close(newvid);
    end
end