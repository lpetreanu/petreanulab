[beads, path] = uigetfile({'W:\shared-petreanu\Camille\CMloop15\20200217\Bleedtrhrough\*.*'},'Load red');
beads = imread([path,'\',beads]);
[gcamp, ~] = uigetfile({[path,'*.*']},'Load green');
gcamp = imread([path,'\',gcamp]);

beads = double(beads);
gcamp = double(gcamp);

beads_2 = reshape(beads, 512*512,1);
gcamp_2 = reshape(gcamp, 512*512,1);

minPixelValue = min(min(min(beads_2)),min(min(gcamp_2)));
beads_2 = beads_2 - minPixelValue;
gcamp_2 = gcamp_2 - minPixelValue;

corr_factor = beads_2\gcamp_2;

figure, hold on
plot(beads_2,gcamp_2,'ko')
xlabel('beads')
ylabel('gcamp')
plot(beads_2, corr_factor*beads_2,'r')
title(['y = ', num2str(corr_factor,3), '*x'])
gcamp_corr = gcamp_2 - beads_2./corr_factor;
gcamp_corr = gcamp_corr + minPixelValue;
gcamp_corr = reshape(gcmap_corr,512,512);

corr_factor = gcamp_2\beads_2;
figure, hold on
plot(gcamp_2,beads_2,'ko')
plot(gcamp_2, corr_factor*gcamp_2,'r')
beads_corr = beads_2 - gcamp_2./corr_factor;
beads_corr = beads_corr + minPixelValue;
beads_corr = reshape(beads_corr,512,512);
beads_corr = int16(beads_corr);

% imwrite(beads_corr,[path 'gcamp_corr.tiff'],'tiff');
tt=Tiff('beads_corr','w');
tagstruct.ImageLength = size(beads_corr,1);
tagstruct.ImageWidth = size(beads_corr,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.SampleFormat=Tiff.SampleFormat.Int;
    tt.setTag(tagstruct);
    tt.write(beads_corr);

    tt.writeDirectory()
close(tt)