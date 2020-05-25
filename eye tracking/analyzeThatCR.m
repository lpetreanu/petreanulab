function [ cx, cy ] = analyzeThatCR(params,varargin)
%[ pd cx py ] = analyzeThatPupil( params,disp,fn )
% analyze image to return corneal reflection center x & y position
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

% cx:   corneal reflection center x position; column 1 is the seed
%       position used to generate the starburst rays, column 2 is fitted position
%       for that frame
% cy:   corneal reflection center y position; same setup as cx


dbug        = 0;    %if 1, then shows you thresholded image, sobel filter and crossing points

cx = nan; cy = nan;

if length(varargin) >= 1
    [eyeMoviePath,name,ext] = fileparts(varargin{1});
    eyeMoviefn = [name,ext];
else
    [eyeMoviefn eyeMoviePath]=uigetfile('.avi','Select eye tracking movie file');
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
cxf         = params.crInit(1);
cyf         = params.crInit(2);
pMax        = params.pMax; 
startFrame  = params.startFrame;
sthresh     = params.sthresh;
pthresh     = params.pthresh;
cfthresh    = params.cfthresh;  
mvthresh    = params.mvthresh;  
filtSize    = params.filtSize;


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
H = fspecial('disk',filtSize);
% H = fspecial('disk',3);

nfr = size(img,3);
cx = nan(nfr,2); cy = nan(nfr,2);
disp('Fitting pupil...')
lastfitind = [];
for f = 1:nfr
    if mod(f,50)== 0
        disp(['frame ' num2str(f)])
    end
        
    %enhance contrast
    if enhance
        img_ec = histeq(img(:,:,f));
    else
        img_ec =img(:,:,f);
    end
    
    blurred = imfilter(img_ec,H,'replicate');

    if dbug
        figure(200),subplot(2,2,1)
        imagesc(img(:,:,1));
        colormap gray
        subplot(2,2,2)
        imagesc(blurred); title('Image, filtered')
    end
    %threshold image
    t = prctile(double(blurred(:)),cfthresh); %threshold
    imgThresh = blurred > t;
    while sum(imgThresh(:)) == 0
        t = t-1;
        imgThresh = blurred > t;
    end
    %for first frame, find center of pupil using center of mass, afterwards
    %use center from previous frame
    if f == 1 && ~exist('cxf','var')
        [cIndy, cIndx] = find(imgThresh == 1);
        cxf = round(mean(cIndx));
        cyf = round(mean(cIndy));
        if isnan(cxf)
            disp('No pupil found. Change settings and try again.')
            break
        end
        
    elseif f>1
        lastfitind = find(~isnan(cx(1:f-1)),1,'last'); %find most recent fit 
        if isempty(lastfitind)
            cxf = params.crInit(1);
            cyf = params.crInit(2);
        else
            cxf = round(cx(lastfitind,2)); cyf = round(cy(lastfitind,2));
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
            hold on
            subplot(2,2,3)
            imagesc(imgThresh), title('thresholded image')
            hold on, scatter(cxf,cyf,20,'r','filled')
            colormap gray
            subplot(2,2,1); hold on
            plot(raysX'+cxf,raysY'+cyf), title('Starburst rays')
            drawnow
        end
        xxx = [];yyy = [];

        raysX((raysX+cxf)>=size(img,2) | (raysX+cxf)<=1)=nan;
        raysY((raysY+cyf)>=size(img,1) | (raysY+cyf)<=1)=nan;
        %finds edges in image
        e = edge(blurred,'sobel',sthresh);
%         ee = imsharpen(double(e));
%         e = ee~=0;
        if dbug
            figure(200), subplot(2,2,4)
            imagesc(e), title('Sobel filter')
        end
        
        [ex,ey] = find(e);
        eInd = sub2ind(size(imgThresh),ex,ey);
        e = double(e);
        
        for i=1:length(angles)
            %             g = zeros(size(imgThresh));
            %             g(raysY(i,:)+pyf,raysX(i,:)+cxf) = 1;
            %             imagesc(g)
            %             [gx,gy] = find(g);
            %             ind = sub2ind(size(imgThresh),gx,gy);
            ind = unique(sub2ind(size(imgThresh),raysY(i,:)+cyf,raysX(i,:)+cxf)); %in positions on the image, in pixels
%             e = double(e);
%             e(ind) = 1;
%             e(eInd) = 2;
%             figure,imagesc(e)
            
            %finds intersection between rays and edges
            intInd = intersect(eInd, ind);
            
            %finds first intersection closest to pupil center
            if ~isempty(intInd)
                [yy,xx]=ind2sub(size(imgThresh),intInd);
                [rayCross,rayXInd]=min(sqrt((xx-cxf).^2+(yy-cyf).^2)); %rayCross is the distance to rayCross in pixels
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
            if xxx(i)-3 < 1 || xxx(i)+3 > size(blurred,2) || yyy(i)-3 < 1 || yyy(i)+3 > size(blurred,1)
                avgLum(i) = pthresh-10;
            else
                avgLum(i) = mean(mean(blurred(yyy(i)-3:yyy(i)+3,xxx(i)-3:xxx(i)+3)));
            end
        end
        
        %kick out outliers
        if dbug
             figure(203);
            eh = imagesc(img_ec); hold on
            colormap gray
            scatter(xxx,yyy,20,'c','filled')
        end
        xxx(avgLum < pthresh) = [];
        yyy(avgLum < pthresh) = [];
        if dbug
            scatter(xxx,yyy,20,'g')
        end
        d = sqrt((xxx-cxf).^2+(yyy-cyf).^2);
        xxx(d>(mean(d)+2*std(d)) | d<(mean(d)-2*std(d))) = [];
        yyy(d>(mean(d)+2*std(d)) | d<(mean(d)-2*std(d))) = [];
        if dbug
            scatter(xxx,yyy,20,'m')
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
                    fitted(f).cx = xc; fitted(f).py = yc; fitted(f).rad = R;
                    fitFound = 1;
                end
            end
        end
        
        %create return variables
        cx(f,1) = cxf; cy(f,1) = cyf;
        if fitFound
            if ellipse
                cd(f) = fitted(f).long_axis;
                cx(f,2) = fitted(f).Y0_in; cy(f,2) = fitted(f).X0_in;
            else
                cd(f) = R*2; cx(f,2) = xc; cy(f,2) = yc;
                %if any parameters change more than mvthresh, throw out the
                %fit        
                if ~isempty(lastfitind)
                    if f > 1 && ((abs(cd(f)-cd(lastfitind)) > mvthresh) || (abs(cx(f,2)-cxf) > mvthresh) || (abs(cy(f,2)-cyf) > mvthresh))
                        cd(f) = nan; cx(f,2) = nan; cy(f,2) = nan;
                    end
                else
                    if f > 1 && ((abs(cx(f,2)-cxf) > mvthresh) || (abs(cy(f,2)-cyf) > mvthresh))
                        cd(f) = nan; cx(f,2) = nan; cy(f,2) = nan;
                    end
                end
            end
        else
            cd(f) = nan;
            cx(f,2) = nan; cy(f,2) = nan;
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
               imagesc(histeq(img(:,:,f)))
           else
               toplot = repmat(histeq(img(:,:,f)),[1,1,3]);
           end
       else
           if slowMovie
               imagesc(img(:,:,f))
           else
               toplot = repmat(img(:,:,f),[1,1,3]);
           end
       end
       if dis == 1 || slowMovie
           title(num2str(f))
           colormap gray
           hold on
           %plots center of pupil
           scatter(cx(f,1),cy(f,1),20,'r','filled')
           scatter(cx(f,2),cy(f,2),25,'b','filled')
       end
        if ~isnan(cd(f)) && cd(f)~= 0
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
                    plot(cos(angles)*cd(f)/2+cx(f,2),sin(angles)*cd(f)/2+cy(f,2))
                    frame = im2frame(getimage(h));
                else
                    
                    % plot circle points
                    try
                        ind = sub2ind(size(toplot),round(sin(angles)*cd(f)/2+cy(f,2)),round(cos(angles)*cd(f)/2+cx(f,2)),3*ones(size(angles)));
                    catch
                        ind = [];
                        disp(['fitted circle out of image range for frame ' num2str(f)])
                    end
                    toplot(ind) = 255;

                    % plot center
                    try
                        toplot(round(cy(f,2))-3:round(cy(f,2))+3,round(cx(f,2))-3:round(cx(f,2))+3,1:3) = 255;
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

