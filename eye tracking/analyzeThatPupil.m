function [ pd, px, py ] = analyzeThatPupil(params,varargin)
%[ pd px py ] = analyzeThatPupil( params,disp,fn )
% analyze image to return pupil diameter, and center x & y position
% inspired by Zoccolan et al 2010 Frontiers in Neuroscience

% INPUTS

% params:   structure with parameters for analysis
%   eyeInit:    initial center of pupil in the format [x y]
%   enhance:    set to 1 to boost contrast; otherwise, leave at 0 or blank
%   s:          fit circle 'c' or ellipse 'e'
%   pMax        max radius of pupil(in pixels)
%   startFrame
%   sthresh     threshold for sobel filter; higher means more restrictive (more restrictive means it finds less edges)
%   pthresh     percent threshold for darkest pixels (to count as pupil)
%   cfthresh    threshold for points to be counted as edge of corneal reflection
%   mvthresh    max movement in pixels allowed from frame to frame
% fn:       name of file to analyze, if you don't want it to prompt
% disp:     0(default) : don't display or save resulted fit as movie
%           1 : display and save movie with fitted pupil
%           2 : save, but don't display movie with fitted pupil
%           movie saved as'fit_ + filename'

% OUTPUTS
% pd:   pupil diameter
% px:   pupil center x position; each row is a frame, column 1 is the seed
%       position used to generate the starburst rays, column 2 is fitted position
%       for that frame
% py:   pupil center y position; same setup as px



dbug        = 0;    %if 1, then shows you thresholded image, sobel filter and crossing points

pd = nan; px = nan; py = nan;

if length(varargin) >= 1
    [eyeMoviePath,name,ext] = fileparts(varargin{1});
    eyeMoviefn = [name,ext];
else
    [eyeMoviefn, eyeMoviePath]=uigetfile('.avi','Select eye tracking movie file');
end

if length(varargin) >= 2
    dis = varargin{2};
else
    dis = 0;
end

%set some parameters
if strcmp(params.s,'e')
    ellipse = 1; circle = 0;
elseif strcmp(params.s,'c')
    ellipse = 0; circle = 1;
else
    disp('Please tell me if I should fit a circle or an ellipse')
    return;
end
enhance     = params.enhance;
pxf         = params.eyeInit(1);
pyf         = params.eyeInit(2);
pMax        = params.pMax;
startFrame  = params.startFrame;
sthresh     = params.sthresh;
pthresh     = params.pthresh;
cfthresh    = params.cfthresh;
mvthresh    = params.mvthresh;
filtSize    = params.filtSize;
if isfield(params,'deltapdthresh')
    deltapdthresh = params.deltapdthresh;
else
    deltapdthresh = mvthresh;
end


slowMovie   = 0;    %makes movies nicer, but slower; NEEDS TO BE 1 FOR FITTING ELLIPSES!
if ellipse; slowMovie = 1; end
if dis == 1; slowMovie = 1; end

%read in video file
vidObj = VideoReader([eyeMoviePath filesep eyeMoviefn]);

endFrame    = vidObj.NumberOfFrames;

%read frames
img = read(vidObj,[startFrame endFrame]);
img = img(:,:,1,:); %take only red channel
img = squeeze(img);

%create filter to smooth out edges
% H = fspecial('disk',filtSize);

nfr = size(img,3);
px = nan(nfr,2); py = nan(nfr,2); pd = nan(nfr,1);
disp('Fitting pupil...')
lastfitind = [];

% initiate figures
if dbug
    f1 = figure('Position',[175 560 560 410],'Name','Check params');
    f2 = figure('Position',[175 60 660 410],'Name','Cyan: edge cross; green: passed CF thresh; magenta: distance thresh');
    f3 = figure('Position',[800 560 900 260],'Name','Check thresholds');
end
for f = 1:nfr
    if mod(f,20)== 0
        disp(['frame ' num2str(f)])
    end
    
    %enhance contrast
    if enhance
        img_ec = histeq(img(:,:,f));
    else
        img_ec =img(:,:,f);
    end
    
    % filter image
    if filtSize>0
        %     blurred = imfilter(img_ec,H,'replicate');
        blurred = imgaussfilt(img_ec,filtSize); % maybe gaussian low pass filter works better?
    else
        blurred = img_ec;
    end
    
    if dbug
        figure(f1)
        subplot(2,2,1)
        imagesc(img(:,:,1)); axis equal
        colormap gray
        subplot(2,2,2)
        imagesc(blurred); title('Image, filtered'); axis equal
        figure(f3)
        subplot(1,2,1)
        imagesc(discretize(double(blurred),10)); axis equal; axis off
        colormap gray; colorbar; title('Discretized into 10 bins (for crthresh)')
        subplot(1,2,2)
        edges = [0 quantile(double(blurred(:)),11) 255];
        imagesc(discretize(double(blurred),edges))
        axis equal; axis off
        colormap gray; colorbar; title('Quantiles (for pthresh)')
    end
    %threshold image
    t = prctile(double(blurred(:)),pthresh); %threshold
    imgThresh = blurred < t;
    while sum(imgThresh(:)) == 0
        t = t+1;
        imgThresh = blurred < t;
    end
    %for first frame, find center of pupil using center of mass, afterwards
    %use center from previous frame
    if f == 1 && ~exist('pxf','var')
        [pIndy, pIndx] = find(imgThresh == 1);
        pxf = round(mean(pIndx));
        pyf = round(mean(pIndy));
        if isnan(pxf)
            disp('No pupil found. Change settings and try again.')
            break
        end
        
    elseif f>1
        lastfitind = find(~isnan(pd(1:f-1)),1,'last'); %find most recent fit
        if isempty(lastfitind)
            pxf = params.eyeInit(1);
            pyf = params.eyeInit(2);
        else
            pxf = round(px(lastfitind,2)); pyf = round(py(lastfitind,2));
        end
    end
    
    
    %starburst algorithm
    %make rays
    angles = 0:.05:2*pi;
    rayTipx = cos(angles)*pMax;
    rayTipy = sin(angles)*pMax;
    for i =1:length(angles)
        raysX(i,:) = round(linspace(0,rayTipx(i),500));%in pixels
        raysY(i,:) = round(linspace(0,rayTipy(i),500));%in pixels
    end
    
    if dbug
        figure(f1)
        hold on
        subplot(2,2,3)
        imagesc(imgThresh); axis equal;title('thresholded image')
        hold on, scatter(pxf,pyf,20,'r','filled')
        colormap gray
        subplot(2,2,1); hold on
        plot(raysX'+pxf,raysY'+pyf), title('Starburst rays')
        drawnow
        waitforbuttonpress()
    end
    xxx = [];yyy = [];
    
    raysX((raysX+pxf)>=size(img,2) | (raysX+pxf)<=1)=nan;
    raysY((raysY+pyf)>=size(img,1) | (raysY+pyf)<=1)=nan;
    %finds edges in image
    e = edge(blurred,'sobel',sthresh);
    %         ee = imsharpen(double(e));
    %         e = ee~=0;
    if dbug
        figure(f1), subplot(2,2,4)
        imagesc(e), axis equal; title('Sobel filter')
    end
    
    [ex,ey] = find(e);
    eInd = sub2ind(size(imgThresh),ex,ey);
    e = double(e);
    
    for i=1:length(angles)
        %             g = zeros(size(imgThresh));
        %             g(raysY(i,:)+pyf,raysX(i,:)+pxf) = 1;
        %             imagesc(g)
        %             [gx,gy] = find(g);
        %             ind = sub2ind(size(imgThresh),gx,gy);
        ind = unique(sub2ind(size(imgThresh),raysY(i,:)+pyf,raysX(i,:)+pxf)); %in positions on the image, in pixels
        %             e = double(e);
        %             e(ind) = 1;
        %             e(eInd) = 2;
        %             figure,imagesc(e)
        
        %finds intersection between rays and edges
        intInd = intersect(eInd, ind);
        
        %finds first intersection closest to pupil center
        if ~isempty(intInd)
            [yy,xx]=ind2sub(size(imgThresh),intInd);
            [rayCross,rayXInd]=min(sqrt((xx-pxf).^2+(yy-pyf).^2)); %rayCross is the distance to rayCross in pixels
            rayCrossX = xx(rayXInd);
            rayCrossY = yy(rayXInd);
            xxx = [xxx rayCrossX];
            yyy = [yyy rayCrossY];
        end
    end
    %get unique values
    intIndU = unique(sub2ind(size(imgThresh),yyy,xxx));
    [yyy,xxx] = ind2sub(size(imgThresh),intIndU);
    
    %check that these are points on the pupil border (not border of
    %reflection)
    avgLum = zeros(size(xxx));
    for i =1:length(xxx)
        %checks if pixel at the edge of image
        if xxx(i)-3 < 1 || xxx(i)+3 > size(blurred,2) || yyy(i)-3 < 1 || yyy(i)+3 > size(blurred,1)
            avgLum(i) = cfthresh+10;
        else
            avgLum(i) = mean(mean(img_ec(yyy(i)-3:yyy(i)+3,xxx(i)-3:xxx(i)+3)));
        end
    end
    
    %kick out outliers
    if dbug
        figure(f2)
        eh = imagesc(img_ec); hold on
        colormap gray; axis equal
        scatter(xxx,yyy,20,'c','filled')
    end
    xxx(avgLum > cfthresh) = [];
    yyy(avgLum > cfthresh) = [];
    if dbug
        figure(f2)
        scatter(xxx,yyy,20,'g')
    end
    d = sqrt((xxx-pxf).^2+(yyy-pyf).^2);
    xxx(d>(mean(d)+2*std(d)) | d<(mean(d)-2*std(d))) = [];
    yyy(d>(mean(d)+2*std(d)) | d<(mean(d)-2*std(d))) = [];
    if dbug
        figure(f2)
        scatter(xxx,yyy,20,'m');
        hold on; scatter(pxf,pyf,20,'r','filled')
        title(sprintf('frame %i',f))
        drawnow
        pause(.3)
    end
    
    if ellipse
        %fit an ellipse
        fitFound = 1;
        lastwarn('')
        if length(xxx)<5
            fitFound = 0;
        else
            if dbug
                fitted(f) = fit_ellipse_m( xxx,yyy, eh );
            else
                fitted(f) = fit_ellipse_m( xxx,yyy);
            end
            if (strcmp(lastwarn,'fit_ellipse: Did not locate an ellipse')) || (strcmp(lastwarn,'Matrix is singular to working precision.'))
                fitFound = 0; %means ellipse not fit
                disp(f)
            end
        end
    elseif circle
        lastwarn('')
        if length(xxx)<5
            fitFound = 0;
        else
            [xc,yc,R] = circfit(xxx,yyy);
            w = warning('query','last');
            if strcmp(lastwarn,'')==0 %strcmp(w.identifier,'MATLAB:rankDeficientMatrix')
                fitFound = 0;
            else
                fitted(f).px = xc; fitted(f).py = yc; fitted(f).rad = R;
                fitFound = 1;
            end
        end
    end
    
    %create return variables
    px(f,1) = pxf; py(f,1) = pyf;
    if fitFound
        if ellipse
            pd(f) = fitted(f).long_axis;
            px(f,2) = fitted(f).Y0_in; py(f,2) = fitted(f).X0_in;
        else
            pd(f) = R*2; px(f,2) = xc; py(f,2) = yc;
            %if any parameters change more than mvthresh, throw out the
            %fit  (pd or distance to new center)
            if ~isempty(lastfitind)
                if f > 1 && ((abs(pd(f)-pd(lastfitind)) > deltapdthresh) || sqrt((px(f,2)-pxf).^2+(py(f,2)-pyf).^2) > mvthresh)
                    pd(f) = nan; px(f,2) = nan; py(f,2) = nan;
                end
            else
                if f > 1 && sqrt((px(f,2)-pxf).^2+(py(f,2)-pyf).^2) > mvthresh
                    pd(f) = nan; px(f,2) = nan; py(f,2) = nan;
                end
            end
        end
    else
        pd(f) = nan;
        px(f,2) = nan; py(f,2) = nan;
    end
    
end %end loop through frames

%display fit and/or save movie
if dis > 0
    disp('Building movie with fit...')
    %play movie with center and diameter printed on video
    if dis == 2
        h = figure('visible','off');
    elseif dis == 1
        h = figure(201);
    end
    
    writerObj = VideoWriter(['fit_',eyeMoviefn]);
    open(writerObj);
    for f = 1:nfr
        %        disp(f)
        if enhance
            if slowMovie
                imagesc(histeq(img(:,:,f))); axis equal; axis off
            else
                toplot = repmat(histeq(img(:,:,f)),[1,1,3]);
            end
        else
            if slowMovie
                imagesc(img(:,:,f)); axis equal; axis off
            else
                toplot = repmat(img(:,:,f),[1,1,3]);
            end
        end
        if dis == 1 || slowMovie
            title(num2str(f))
            colormap gray
            hold on
            %plots center of pupil
            scatter(px(f,1),py(f,1),20,'r','filled')
            scatter(px(f,2),py(f,2),25,'b','filled')
        end
        if ~isnan(pd(f)) && pd(f)~= 0
            if ellipse
                R = [cos(fitted(f).phi) sin(fitted(f).phi); -sin(fitted(f).phi) cos(fitted(f).phi) ];
                
                % the axes
                ver_line        = [ [fitted(f).X0 fitted(f).X0]; fitted(f).Y0+fitted(f).b*[-1 1] ];
                horz_line       = [ fitted(f).X0+fitted(f).a*[-1 1]; [fitted(f).Y0 fitted(f).Y0] ];
                new_ver_line    = R*ver_line;
                new_horz_line   = R*horz_line;
                
                % the ellipse
                theta_r         = linspace(0,2*pi);
                ellipse_x_r     = fitted(f).X0 + fitted(f).a*cos( theta_r );
                ellipse_y_r     = fitted(f).Y0 + fitted(f).b*sin( theta_r );
                rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
                % draw
                hold_state = get( gca,'NextPlot' );
                set( gca,'NextPlot','add' );
                plot( new_ver_line(2,:),new_ver_line(1,:),'r' );
                plot( new_horz_line(2,:),new_horz_line(1,:),'r' );
                plot( rotated_ellipse(2,:),rotated_ellipse(1,:),'r' );
                set( gca,'NextPlot',hold_state );
                frame = im2frame(getCdata(h));
            else
                if slowMovie
                    plot(cos(angles)*pd(f)/2+px(f,2),sin(angles)*pd(f)/2+py(f,2),'LineWidth',2)
                    %                     frame = im2frame(getimage(h));
                    frame = getframe(h);
                else
                    
                    % plot circle points
                    try
                        ind = sub2ind(size(toplot),round(sin(angles)*pd(f)/2+py(f,2)),round(cos(angles)*pd(f)/2+px(f,2)),3*ones(size(angles)));
                    catch
                        ind = [];
                        disp(['fitted circle out of image range for frame ' num2str(f)])
                    end
                    toplot(ind) = 255;
                    
                    % plot center
                    try
                        toplot(round(py(f,2))-3:round(py(f,2))+3,round(px(f,2))-3:round(px(f,2))+3,1:3) = 255;
                    catch
                        disp(['Center of fitted pupil is out of range'])
                    end
                    frame = im2frame(toplot);
                end
            end
            
            %drawnow
            try
                writeVideo(writerObj,frame);
            catch ME
                if (strcmp(ME.identifier,'MATLAB:audiovideo:VideoWriter:invalidDimensions'))
                    size(toplot)
                    close(writerObj)
                end
            end
            if dis == 1
                pause(.3)
            end
        end
    end
    close(writerObj);
    set(0,'DefaultFigureVisible','on');
end

%function end
end

