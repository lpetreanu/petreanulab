saveFig_Flag = 1;
saveDir = [data.Info.mainDir, '\Figures\BeadsCells'];
if ~exist(saveDir)
    mkdir(saveDir)
end

Vis_Laser_mean = mean(data.visNopto.Vis_Laser,4);
Vis_Laser_sem = std(data.visNopto.Vis_Laser,[],4)/sqrt(size(data.visNopto.Vis_Laser,4));

Vis_NoLaser_mean = mean(data.visNopto.Vis_NoLaser,4);
Vis_NoLaser_sem = std(data.visNopto.Vis_NoLaser,[],4)/sqrt(size(data.visNopto.Vis_NoLaser,4));

timeVect = data.Info.timeVect;
xvect = data.Info.xvect;

for i  = 1:size(data.VisResponsiveROIs.beads,2)
    if isnan(data.VisResponsiveROIs.beads(2,i))
    else
    ymax = max(max(max(Vis_Laser_mean(data.VisResponsiveROIs.beads(1,i),:,:)+Vis_Laser_sem(data.VisResponsiveROIs.beads(1,i),:,:))),max(max(Vis_NoLaser_mean(data.VisResponsiveROIs.beads(1,i),:,:)+Vis_NoLaser_sem(data.VisResponsiveROIs.beads(1,i),:,:))));
    ymin = min(min(min(Vis_Laser_mean(data.VisResponsiveROIs.beads(1,i),:,:)-Vis_Laser_sem(data.VisResponsiveROIs.beads(1,i),:,:))),min(min(Vis_NoLaser_mean(data.VisResponsiveROIs.beads(1,i),:,:)-Vis_NoLaser_sem(data.VisResponsiveROIs.beads(1,i),:,:))));
    figure;
    for ii = 1:data.Info.nVisStimTypes
        subplot(2,data.Info.nVisStimTypes/2,ii); hold on
        if data.VisResponsiveROIs.beads(2,i) == 1
            h1 = fill([timeVect(xvect) flip(timeVect(xvect))],...
                [Vis_Laser_mean(data.VisResponsiveROIs.beads(1,i),xvect,ii)+Vis_Laser_sem(data.VisResponsiveROIs.beads(1,i),xvect,ii) flip(Vis_Laser_mean(data.VisResponsiveROIs.beads(1,i),xvect,ii)-Vis_Laser_sem(data.VisResponsiveROIs.beads(1,i),xvect,ii))],[1 0 0],'LineStyle','none');
            plot(timeVect(xvect),Vis_Laser_mean(data.VisResponsiveROIs.beads(1,i),xvect,ii),'r')
        elseif data.VisResponsiveROIs.beads(2,i) == 0
            h1 = fill([timeVect(xvect) flip(timeVect(xvect))],...
                [Vis_Laser_mean(data.VisResponsiveROIs.beads(1,i),xvect,ii)+Vis_Laser_sem(data.VisResponsiveROIs.beads(1,i),xvect,ii) flip(Vis_Laser_mean(data.VisResponsiveROIs.beads(1,i),xvect,ii)-Vis_Laser_sem(data.VisResponsiveROIs.beads(1,i),xvect,ii))],[0 0 1],'LineStyle','none');
            plot(timeVect(xvect),Vis_Laser_mean(data.VisResponsiveROIs.beads(1,i),xvect,ii),'b')
        end
        set(h1,'facealpha',.5)
        h2 = fill([timeVect(xvect) flip(timeVect(xvect))],...
            [Vis_NoLaser_mean(data.VisResponsiveROIs.beads(1,i),xvect,ii)+Vis_NoLaser_sem(data.VisResponsiveROIs.beads(1,i),xvect,ii) flip(Vis_NoLaser_mean(data.VisResponsiveROIs.beads(1,i),xvect,ii)-Vis_NoLaser_sem(data.VisResponsiveROIs.beads(1,i),xvect,ii))],[0 0 0],'LineStyle','none');
        set(h2,'facealpha',.5)
        plot(timeVect(xvect),Vis_NoLaser_mean(data.VisResponsiveROIs.beads(1,i),xvect,ii),'k')
        ylim([ymin-0.05 ymax+0.05])
        
        plot([data.Info.xstim_opto],[ymax ymax],'r','LineWidth',5)
        plot([data.Info.xstim_vis],[ymax+0.05 ymax+0.05],'k','LineWidth',5)
        %             plot([data.Info.xstim_opto(1) data.Info.xstim_opto(1)],[ymin ymax],'k:')
        %             plot([data.Info.xstim_opto(2) data.Info.xstim_opto(2)],[ymin ymax],'k:')
        
        if ii == data.Info.nVisStimTypes/2+1
            xlabel('Time from photostim onset (s)');
            ylabel('dF/F')
        end
        %                 blah = find(i == row(:));
        %                 if any(col(blah) == ii)
        %                     title(num2str(p_all(i,ii)));
        %                     box on
        %                 end
        
    end
    sgtitle(['ROI#', num2str(data.VisResponsiveROIs.beads(1,i))]);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.04, 0.80, 0.5]);
    if saveFig_Flag
                saveas(gcf,[saveDir,'\',num2str(data.VisResponsiveROIs.beads(1,i)),'_Directions_Tuning.tif'])
                pause(0.2); close
    end
    end
end