% script to improve contrasts of movies before fitting in DeepLabCut
% Marina Feb 2021
%% set directory and get file names
eyeMoviePath= "E:\Imaging2AFCPop\eye tracking\MF183";
cd(eyeMoviePath)
fn = dir('MF*.avi');
nMov =length(fn);

%% loop through movies
for m = 1:nMov
    eyeMoviefn=fn(m).name;
    disp(eyeMoviefn)
    
    vr = VideoReader(fullfile(eyeMoviePath,eyeMoviefn));
    
    % read frames
    vid = read(vr);
    
    % get number of frames
    [~,~,~,nfr] = size(vid);
    
    %initialize new movie structure
    new = struct('cdata',[],'colormap',[]);
    
    % loop through frames
    for f = 1:nfr
        if mod(f,100)== 0
            disp(['frame ' num2str(f)])
        end
        
        % convert image to grayscale
        gray = rgb2gray(vid(:,:,:,f));
        
        %enhance contrast
        ec = histeq(gray);
        
        %convert back to RGB and save frame
        new(f) = im2frame(repmat(ec,[1,1,3]));
    end
    
    % write new movie
    vw = VideoWriter(strcat('ec_', eyeMoviefn));
    vw.FrameRate = vr.FrameRate;
    open(vw);
    writeVideo(vw,new)
    close(vw)
end