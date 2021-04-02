%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveFig = true; 
selection = 1; % 1=beads-based
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n')
if saveFig
    mainDir =  data.Info.mainDir;
    if exist(mainDir,'dir')
        disp(mainDir)
    else
        mainDir = uigetdir('X:\camille.mazo\2P_processed','Remap main dir');
        fprintf('Main directory has changed to %s \n',mainDir)
        data.Info.mainDir = mainDir;
    end
    saveDir = [data.Info.mainDir, '\Figures\PopPlots'];
    if ~ exist(saveDir)
        mkdir(saveDir)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timeStimInd = data.Info.timeStimInd;
meanResp = squeeze(nanmean(nanmean(data.data_sorted(:,timeStimInd,2:9,:,:,:),6),2));
figure; hold on
for i = 1:size(data.data_sorted,1)
    [best_RespSF1, I1] = max(squeeze(meanResp(i,:,1,1)),[],'omitnan');
    [best_RespSF2, I2] = max(squeeze(meanResp(i,:,2,1)),[],'omitnan');
    if max(best_RespSF1,best_RespSF2) == best_RespSF1
        R = I1;
        C = 1;
    else
        R = I2;
        C = 2;
    end
    plot(meanResp(i,R,C,1),meanResp(i,R,C,2),'ko')
    BestVisStim(i) = meanResp(i,R,C,1);
    BestVisStim_Photostim(i) = meanResp(i,R,C,2);
end
xl = xlim; yl = ylim;
xmax = max(abs(xl(1)),xl(2)); ymax = max(abs(yl(1)),yl(2));
axismax = max(xmax,ymax);
plot([-axismax axismax],[-axismax  axismax ],'k:');
clear xlim ylim
xlim(xl); ylim(yl)
xlabel({'dF/F','Best Vis Stim'}); ylabel({'Same Vis Stim+Photostim','dF/F'})
saveas(gcf,[saveDir '\ImpactOfLaserBestVis_Scatter'],'fig');
screen2png([saveDir '\ImpactOfLaserBestVis_Scatter']);
    
ImpactOfLaser_Vis = BestVisStim_Photostim - BestVisStim;
ImpactOfLaser_Vis_mean = mean(ImpactOfLaser_Vis);
ImpactOfLaser_Vis_sem = std(ImpactOfLaser_Vis)./sqrt(size(ImpactOfLaser_Vis,2));
for i = 1:size(data.data_sorted,1)
    if ImpactOfLaser_Vis(1,i)>0
        pos(i) = true;
    elseif ImpactOfLaser_Vis(1,i)<0
         pos(i) = false;
    end
end
ImpactOfLaser_Vis_Pos_mean = mean(ImpactOfLaser_Vis(pos));
ImpactOfLaser_Vis_Neg_mean = mean(ImpactOfLaser_Vis(~pos));
ImpactOfLaser_Vis_Pos_sem = std(ImpactOfLaser_Vis(pos))./sqrt(size(ImpactOfLaser_Vis(pos),2));
ImpactOfLaser_Vis_Neg_sem = std(ImpactOfLaser_Vis(~pos))./sqrt(size(ImpactOfLaser_Vis(~pos),2));
figure; hold on
bar([1 2 3],[ImpactOfLaser_Vis_mean ImpactOfLaser_Vis_Pos_mean ImpactOfLaser_Vis_Neg_mean])
plot([1 1],[ImpactOfLaser_Vis_mean-ImpactOfLaser_Vis_sem ImpactOfLaser_Vis_mean+ImpactOfLaser_Vis_sem],'k-')
plot([2 2],[ImpactOfLaser_Vis_Pos_mean-ImpactOfLaser_Vis_Pos_sem ImpactOfLaser_Vis_Pos_mean+ImpactOfLaser_Vis_Pos_sem],'k-')
plot([3 3],[ImpactOfLaser_Vis_Neg_mean-ImpactOfLaser_Vis_Neg_sem ImpactOfLaser_Vis_Neg_mean+ImpactOfLaser_Vis_Neg_sem],'k-')
xticks([1 2 3]); xticklabels({'All','Pos','Neg'})
title('Impact of Laser on Best Vis Resp')
saveas(gcf,[saveDir '\ImpactOfLaserBestVis_Bar'],'fig');
screen2png([saveDir '\ImpactOfLaserBestVis_Bar']);

meanResp_spont = squeeze(nanmean(nanmean(data.data_sorted(:,timeStimInd,1,1,:,:),6),2));
SpontActivity = squeeze(meanResp_spont(:,1));
PhotostimOnly = squeeze(meanResp_spont(i,2));
figure; hold on
for i = 1:size(data.data_sorted,1)
    plot(squeeze(meanResp_spont(i,1)),squeeze(meanResp_spont(i,2)),'ko');
end
xl = xlim; yl = ylim;
xmax = max(abs(xl(1)),xl(2)); ymax = max(abs(yl(1)),yl(2));
axismax = max(xmax,ymax);
plot([-axismax axismax],[-axismax  axismax ],'k:');
clear xlim ylim
xlim(xl); ylim(yl)
xlabel({'dF/F','Spont Activity'}); ylabel({'Photostim','dF/F'})
saveas(gcf,[saveDir '\ImpactOfLaserSpont_Scatter'],'fig');
screen2png([saveDir '\ImpactOfLaserSpont_Scatter']);

ImpactOfLaser_Spont = PhotostimOnly - SpontActivity;