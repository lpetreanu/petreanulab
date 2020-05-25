function [ap] = pupilPixelsToDegrees(CR,mmPerPix)
%pupilPixelsToDegrees converts eye position and pupil diameter from pixels
%to degrees, based on Stahl et al 2000
%  conversion based on following equation E=arcsin{(CR_P)/Rp }
%   CR is the linear position, P is the video pipil, E is the angular
%   position of the eye, Rp is the radius of the pupil, approximated to
%   be 1mm (based on this paper)
Rp = 1;
ap = asin(CR*mmPerPix/Rp)*180/pi;
end

