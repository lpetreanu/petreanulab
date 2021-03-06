Here is some code for doing rigid registration using only a cropped portion of the image, that is then applied to the whole FOV. 

This is useful for when there is uneven light contamination at the edges of the images. This can happen thanks the to the circuit that shuts down the monitor during scanning so that monitor bleedthrough only happens when the mirrors are turning and flyback. Since this signal is quite strong and does not move with the brain, it can dominate the signal and confuse your registration.

The code was developed for SI4, probably needs a makeover for the latest SI versions. 

You can change the cropping window in Marina_registration_MultPeak_singleTrial_padding_multPlanes.m

Example (frame-by-frame registration):
`Marina_registration_MultPeaks_acrossTrials_MultPlanes('Reg_',1,[1 2 3 4 5], 5, 'SI4',1)`

Can also do registration on a trial-by-trial basis (more suited to low SNR, anesthetized recordings where movement is expected to primarily be due to drift),
just change the last argument (byFrame) to 0. 

### To create registration target:
**trial average**
ImageJ: 
- concatenate several trials (Image >> Stacks >> Tools >> Concatenate)
- run extractTargetMultPlane_MF.txt

**registered trial average**

``` load('Reg_MF379_240818_Loc2_MaxProjection_Plane_1.mat')
p1 = mean(MaxProjectionPlane,3);
figure,imagesc(p1)
colormap gray; cmap = colormap;
imwrite(p1,cmap,'Pl1.tif')
``` 

**registered trial average with subset of trials**
```
p1 = mean(MaxProjectionPlane(:,:,[4,6,8,9]),3);
```
