% Prepare eye videos from face recordings in quietFB experiment 
% 1) fix file name to reflect trial # (+1) and format trial #
% 2) crop file to only include eye
% 3) save file 
clearvars
%% set group to analyze
fr = 8.5;
co = [0 0.35 0.57;
    0.12 0.63 0.67]; % color scheme
%set some things
prot = 'natImgs';
plt = 0; %whether to plot
%% session files (missing GF337 20190619/20 bc no dark protocol that day)
dr = 'E:\quietFB\';
cd(dr)
mice = {'GF337','GF337','GF338','GF338','GF339','GF339','GF339',...
    'GF341','GF343','GF343','GF343','GF344','GF344','GF346','GF346','GF346',...
    'GF320','GF321','GF322','GF330','GF335','GF334','GF335','GF335','GF335'};
cno_day = {'20190618','20190620','20190618','20190702',...
    '20190716','20190802','20190730','20190828','20190903','20190912',...
    '20190917','20190903','20190912','20190927','20191001','20191004',...
    '20190201','20190201','20190131','20190426','20190510','20190524',...
    '20190514','20190521','20190524'};
sal_day = {'20190617',...
    '20190619',...
    '20190617',...
    '20190701',...
    '20190715','20190801','20190729','20190827','20190902','20190911',...
    '20190916','20190902','20190911','20190926','20190930','20191003',...
    '20190130','20190129','20190129','20190425','20190509','20190523',...
    '20190513','20190520','20190523'};

%% loop through all sessions

for m = 1:length(mice)
    mouse = mice{m};
    %process saline day of recordings
    dt = sal_day{m};
    doMagic(mouse, dt, prot)
%    
    %process CNO day of recordings
    dt = cno_day{m};
    doMagic(mouse, dt, prot)
end
%%
function doMagic(mouse, dt, prot)
% set directory and base name
d = sprintf('E:\\quietFB\\faceVids\\%s\\%s\\%s',mouse,circshift(dt,4),prot);
base = sprintf('%s_%s_%s',mouse,circshift(dt,4),prot);

cd(d)
fns = dir('*.avi');
% makes new directory where to put cropped eye videos
new_d = sprintf('E:\\quietFB\\eyeVids\\%s\\%s\\%s',mouse,circshift(dt,4),prot);
mkdir(new_d)

% loops through movie files
for t = 1:length(fns)
    % make old and new filenames
    oldFn = sprintf('%s%d.avi',base,t-1);
    newFn = sprintf('%s_%03d.avi',base,t);
    fprintf('%s...',oldFn)
    
    % open face videofile
    faceVid = VideoReader([d filesep oldFn]);
    
    % read all frames
    f=1; allFrames = [];
    while hasFrame(faceVid)
        allFrames(:,:,:,f) = readFrame(faceVid);
        f = f+1;
    end
    allFrames = allFrames(:,:,1,:); %take only red channel
    
    % find crop box (1st trial only)
    if t == 1
        figure,imagesc(squeeze(allFrames(:,:,1,1))); colormap gray; axis equal; axis off
        roiInd = round(getrect(gcf)); %returns indices of roi to analyze, [xmin ymin width height]
        str = input('Are you satisfied with the cropbox you drew? [y/n]','s');
        while str=='n'
            roiInd = round(getrect(gcf)); %returns indices of roi to analyze, [xmin ymin width height]
            str = input('Are you satisfied with the cropbox you drew? [y/n]','s');
        end
        close(gcf)
    end
    
    %crop video
    roi = allFrames(roiInd(2):roiInd(2)+roiInd(4),roiInd(1):roiInd(1)+roiInd(3),:,:);
    
    % rescale pixel values to [0,1] so that VideoWriter doesn't complain
    [a,b,c,e] = size(roi);
    roi = reshape(rescale(roi(:)),[a,b,c,e]);

    % save cropped video
    eyeVid = VideoWriter([new_d filesep newFn]);
    open(eyeVid);
    writeVideo(eyeVid,roi)
    close(eyeVid);
    fprintf('done\n')
end
end