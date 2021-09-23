function im_smth = im_mov_avg(im, span)

% Moving average of images over a series of frams.
% Note, span has to be odd number

% - NX 3/11/2009 LTP2012

mean_im = mean(im,3);

pad=repmat(mean_im,[1 1 (span-1)/2]);

temp = cat(3,pad,im,pad);

im_smth = zeros(size(im));

for i = 1:size(im,3) % (span-1)/2+1 : size(im,3)+(span-1)/2
    im_smth(:,:,i) = mean(temp(:,:,i:i+span-1), 3);
end;
im_smth = uint16(im_smth);
